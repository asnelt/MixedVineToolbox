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

function [info,stderr] = mixedvineinfo(vines,pcond,alpha,erreps,cases)
% MIXEDVINEINFO Mixed copula vine mutual information estimate.
%   [INFO,STDERR] = MIXEDVINEINFO(VINES,PCOND,ALPHA,ERREPS,CASES) returns
%   an estimate INFO of the mutual information between a discrete random
%   variable and a set of conditional mixed vines specified in the cell
%   VINES. Each element of this cell specifies as mixed vine as returned by
%   MIXEDVINEFIT. The probability of occurence of each vine is specified in
%   the vector PCOND which must have the same size as VINES.
%   ALPHA is the significance level of the estimate (default ALPHA = 0.05).
%   ERREPS is the maximum standard error of the mutual information estimate
%   as a stopping criterion (default ERREPS = 1e-3).
%   CASES is the number of samples that are drawn in each iteration of the
%   Monte Carlo estimation (default CASES = 1000).
%
%   INFO and STDERR are scalars in unit bit (base 2 logarithm).

% Argument checks
if nargin < 2
    error('mixedvineinfo: Usage [info,stderr] = mixedvineinfo(vines,pcond,alpha,erreps,cases)');
end
if ~iscell(vines)
    error('mixedvineinfo: Argument "vines" must be a cell array');
end
if ~isvector(pcond) || any(pcond<0) || any(pcond>1)
    error('mixedvineinfo: pcond must be a vector of probabilities');
end
if nargin < 3 || isempty(alpha)
    alpha = 0.05;
elseif ~isscalar(alpha)
    error('mixedvineinfo: Argument "alpha" must be a scalar');
end
if nargin < 4 || isempty(erreps)
    erreps = 1e-3;
elseif ~isscalar(erreps)
    error('mixedvineinfo: Argument "erreps" must be a scalar');
end
if nargin < 5
    cases = 1000;
elseif ~isscalar(cases)
    error('mixedvineinfo: Argument "cases" must be a scalar');
end

% Number of conditions
ncond = length(pcond);
% Probabilities of conditions
pcond = repmat(pcond(:)',cases,1);

% Gaussian confidence interval for erreps and level alpha
conf = norminv(1 - alpha,0,1);

% Estimate unconditional entropy h
stderr = inf;
h = 0;
varsum = 0;
k = 0;
x = zeros(cases,length(vines{1}.margins));
pcum = cumsum(pcond(:,end:-1:1),2);
pcum = pcum(:,end:-1:1);
while stderr >= erreps/2
    % Generate samples
    c = repmat(rand(cases,1),1,ncond);
    c = sum(c <= pcum,2);
    for i = 1:ncond
        sel = c==i;
        cases_sel = sum(sel);
        x(sel,:) = mixedvinernd(vines{i},cases_sel);
    end
    p = zeros(cases,1);
    for i = 1:ncond
        p = p + pcond(1,i) * mixedvinepdf(vines{i},x);
    end
    logp = log(p);
    log2p = logp(~isinf(logp)) / log(2);
    k = k + 1;
    % Monte-Carlo estimate of entropy
    h = h + (-mean(log2p) - h) / k;
    % Estimate standard error
    varsum = varsum + sum((-log2p - h) .^ 2);
    stderr = conf * sqrt(varsum / (k * cases * (k * cases - 1)));
end

% Subtract conditional entropies
info = h;
for i = 1:ncond
    errepscond = pcond(1,i) * erreps / 2;
    [hcond,stderrcond] = mixedvineentropy(vines{i},alpha,errepscond,cases);
    info = info - pcond(1,i) * hcond;
    stderr = stderr + pcond(1,i) * stderrcond; 
end

end
