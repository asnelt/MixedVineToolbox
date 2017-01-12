% Copyright (C) 2016 Arno Onken
%
% This file is part of the Mixed Vine Toolbox.
%
% The Mixed Vine Toolbox is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License as published
% by the Free Software Foundation; either version 3 of the License, or (at
% your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
% Public License for more details.
%
% You should have received a copy of the GNU General Public License along
% with this program; if not, see <http://www.gnu.org/licenses/>.

function [vine,logp] = mixedvinefit(x,type,iscont,tl,refine)
% MIXEDVINEFIT Mixed copula vine estimates.
%   VINE = MIXEDVINEFIT(X,TYPE,ISCONT,TL,REFINE) selects the best fitting
%   mixed margin copula vine of type TYPE given the data X, where X has the
%   size [N, D] for N the number of samples and D the dimension of each
%   sample. The selection criterion is the Akaike information criterion for
%   both margins and copula families. TYPE can only be 'c-vine' for the
%   canonical vine for now.
%   The mixed continuous and discrete margins are returned in the cell
%   field VINE.margins. Each element of the cell specifies one margin and
%   is  specified in a struct as returned by MARGINFIT. The margin
%   distribution is continuous if the corresponding element in the boolean
%   vector ISCONT is true and discrete otherwise. TL is a scalar specifing
%   the truncation level of the vine. Elements in the vine tree beyond tree
%   level TL are assumed to be independent. REFINE is a boolean specifying
%   whether the inference for margins results should be refined by means of
%   joint parameter estimation (default REFINE = D<=5).
%   The vine copula families are returned in the D x D cell field
%   VINE.families. Each element of this cell specifies a copula family and
%   can be one of
%    'ind'           for independence,
%    'gaussian'      for the Gaussian copula family,
%    'student'       for the student copula family,
%    'clayton'       for the Clayton copula family,
%    'claytonrot090' for the 90° clockwise rotated Clayton copula family,
%    'claytonrot180' for the survival Clayton copula family,
%    'claytonrot270' for the 270° clockwise rotated Clayton copula family,
%   or be empty if the vine type does not use the element. The parameters
%   of the copula families are returned in the corresponding elements of
%   the cell field VINE.theta.
%   The return value logp is the log likelihood of X given VINE.

% Argument checks
if nargin < 3
    error('mixedvinefit: Usage vine = mixedvinefit(x,type,iscont,tl,refine)');
end
if ~ismatrix(x)
    error('mixedvinefit: Argument "x" must be a matrix');
end
% Check whether the variance in any margin is 0
if any(var(x)==0)
    error('mixedvinefit: Zero variance margin');
end
[cases,d] = size(x);

global newvine newnode u truncation;

if ~ischar(type)
    error('mixedvinefit: Argument "type" must be a string');
end
if ~strcmpi(type,'c-vine')
    error('mixedvinefit: Currently, argument "type" must be "c-vine"');
end
if ~islogical(iscont) || length(iscont) ~= d
    error('mixedvinefit: Argument "iscont" must be a boolean vector of length d');
end
if nargin < 4 || isempty(tl)
    truncation = d-1;
else
    if ~isscalar(tl) || round(tl)~=tl
        error('mixedvinefit: Truncation level "tl" must be an integer scalar');
    end
    truncation = tl;
end
if nargin < 5
    refine = d<=5;
end

newnode = true(d);
u = zeros(cases,d,d);
% Fit margins
margins = cell(d,1);
for i = 1:d
    margins{i} = marginfit(x(:,i),iscont(i));
    u(:,i,i) = margincdf(margins{i},x(:,i));
end
newvine.margins = margins;
newvine.type = type;
newvine.families = cell(d);
newvine.theta = cell(d);

% Fit vine
fitctree(d-1,d);
vine = newvine;
clear global newvine newnode u truncation;

% Joint parameter estimation with current values as initial values
jointtheta = jointpar(vine.theta);
if ~isempty(jointtheta)
    [lb,ub] = vinebounds(vine);
    % Objective function
    objf = @(jointtheta) -sum(log(mixedvinepdf(thetavine(vine,jointtheta),x)+eps));
    % Check whether the initial parameters are feasible
    val = objf(jointtheta);
    while isnan(val) || isinf(val)
        % Try distortion of initpar
        jointtheta = jointtheta + rand(size(jointtheta));
        % Ensure bounds
        k = jointtheta <= lb;
        jointtheta(k) = lb(k) + 1e-3;
        k = jointtheta >= ub;
        jointtheta(k) = ub(k) - 1e-3;
        val = objf(jointtheta);
    end
    if refine
        options = optimset('Algorithm','interior-point','Display','off','MaxIter',100);
        try
            [jointtheta,val] = fmincon(objf,jointtheta,[],[],[],[],lb,ub,[],options);
        catch err
            disp(['mixedvinefit: Unable to refine parameters.' err.message]);
        end
    end
    vine = thetavine(vine,jointtheta);
else
    val = -sum(log(mixedvinepdf(vine,x)+eps));
end
logp = -val;

end


function fitnode(v,i,j,treelevel)
% Fits a single node
global newvine truncation;

if treelevel > truncation
    newvine.families{i,j} = 'ind';
else
    families = {'ind','gaussian','student','clayton','claytonrot090','claytonrot180','claytonrot270'};
    aic = zeros(length(families),1);
    theta = cell(length(families),1);
    for k = 1:length(families)
        theta{k} = copulafit(families{k},v);
        logp = sum(log(copulapdf(families{k},v,theta{k})+eps));
        aic(k) = 2*length(theta{k}) - 2*logp;
    end
    % Use AIC for selecting the best family
    [~,kmin] = min(aic);
    newvine.families{i,j} = families{kmin};
    newvine.theta{i,j} = theta{kmin};
end

end


function fitctree(i,j)
% Fits a canonical vine tree
global newvine newnode u;

if i == 1
    % Leaf
    v = [u(:,1,1) u(:,j,j)];
else
    % Node
    if newnode(i-1,i)
        fitctree(i-1,i);
    end
    if newnode(i-1,j)
        fitctree(i-1,j);
    end
    v = [u(:,i-1,i) u(:,i-1,j)];
end
treelevel = i;
fitnode(v,i,j,treelevel);
u(:,i,j) = copulaccdf(newvine.families{i,j},v,newvine.theta{i,j},1);
newnode(i,j) = false;

end


function jointtheta = jointpar(theta)
% Puts all parameters into a single vector
npar = 0;
for i = 1:size(theta,1)
    for j = 1:size(theta,2)
        npar = npar + length(theta{i,j});
    end
end
jointtheta = zeros(npar,1);
index = 1;
for i = 1:size(theta,1)
    for j = 1:size(theta,2)
        for k = 1:length(theta{i,j})
            jointtheta(index) = theta{i,j}(k);
            index = index + 1;
        end
    end
end

end


function vine = thetavine(vine,jointtheta)
% Sets VINE.theta with the complete parameter vector JOINTTHETA
theta = cell(size(vine.theta));
index = 1;
for i = 1:size(vine.theta,1)
    for j = 1:size(vine.theta,2)
        theta{i,j} = zeros(size(vine.theta{i,j}));
        for k = 1:length(theta{i,j})
            theta{i,j}(k) = jointtheta(index);
            index = index + 1;
        end
    end
end
vine.theta = theta;

end


function [lb,ub] = vinebounds(vine)
% Returns the lower and upper bounds of the vine parameters
npar = 0;
for i = 1:size(vine.theta,1)
    for j = 1:size(vine.theta,2)
        npar = npar + length(vine.theta{i,j});
    end
end
lb = zeros(npar,1);
ub = zeros(npar,1);
index = 1;
for i = 1:size(vine.theta,1)
    for j = 1:size(vine.theta,2)
        if ~isempty(vine.theta{i,j})
            switch lower(vine.families{i,j})
                case 'gaussian'
                    lb(index) = -1;
                    ub(index) = 1;
                case 'student'
                    lb(index) = -1;
                    ub(index) = 1;
                    lb(index+1) = 1e-1;
                    ub(index+1) = 1000;
                case {'clayton','claytonrot090','claytonrot180','claytonrot270'}
                    lb(index) = 1e-3;
                    ub(index) = 20;
            end
            index = index + length(vine.theta{i,j});
        end
    end
end

end
