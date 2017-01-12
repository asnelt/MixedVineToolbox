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

function [p,logp] = mixedvinepdf(vine,x)
% MIXEDVINEPDF Mixed copula vine probability density function.
%   [P,LOGP] = MIXEDVINEPDF(VINE,X) returns the probability density in P
%   and the natural logarithm of the probability density in LOGP of the
%   mixed copula vine specified in the struct VINE at the values in X. X
%   has size [N, D], where N is the number of samples and D is the
%   dimension of each sample. The vine type is specified in the field
%   VINE.type which can only be one 'c-vine' for the canonical vine for
%   now.
%   The mixed continuous and discrete margins are specified in cell field
%   VINE.margins. Each element of the cell specifies one margin and is
%   specified in a struct as returned by MARGINFIT. The vine copula
%   families are specified in the D x D cell field VINE.families. Each
%   element of this cell specifies a copula family and can be one of
%    'ind'           for independence,
%    'gaussian'      for the Gaussian copula family,
%    'student'       for the student copula family,
%    'clayton'       for the Clayton copula family,
%    'claytonrot090' for the 90° clockwise rotated Clayton copula family,
%    'claytonrot180' for the survival Clayton copula family,
%    'claytonrot270' for the 270° clockwise rotated Clayton copula family,
%   or is empty if the vine type does not use the element. The parameters
%   of the copula families are specified in the corresponding elements of
%   the cell field VINE.theta.
%
%   P and LOGP are N-dimensional vectors.

% Argument checks
if nargin < 1
    error('mixedvinepdf: Usage [p,logp] = mixedvinepdf(vine,x)');
end
if ~isstruct(vine)
    error('mixedvinepdf: Argument "vine" must be a struct');
end
d = length(vine.margins);
if ~ismatrix(x)
    error('mixedvinepdf: Argument "x" must be a matrix');
end
if size(x,2) ~= d
    error('mixedvinepdf: Second dimension of u must match number of margins');
end
cases = size(x,1);

global newnode Fp Fm logf;
newnode = true(d);
Fp = zeros(cases,d,d);
Fm = zeros(cases,d,d);
% logf(s,i,j) represents the s-th sample of log(f(j|i))
logf = zeros(cases,d,d);
for i = 1:d
    Fp(:,i,i) = margincdf(vine.margins{i},x(:,i));
    if vine.margins{i}.iscont
        logf(:,i,i) = log(marginpdf(vine.margins{i},x(:,i)));
    else
        Fm(:,i,i) = margincdf(vine.margins{i},x(:,i)-1);
        logf(:,i,i) = log(Fp(:,i,i) - Fm(:,i,i));
    end
end
evalctree(vine,d-1,d);
logp = logf(:,1,1);
for i = 2:d
    logp = logp + logf(:,i-1,i);
end
logp = reshape(logp,[cases 1]);
% Correct numerical inaccuracies
logp = real(logp);
% Massless intervals can result in NaNs; set to 0
logp(isnan(logp)) = -inf;
p = exp(logp);

clear global newnode Fp Fm logf;

end


function evalrightnode(vine,i,j,k,l)

global Fp Fm logf;

if vine.margins{i}.iscont && vine.margins{j}.iscont
    Fp(:,i,j) = copulaccdf(vine.families{j,i},[Fp(:,l,j) Fp(:,k,i)],vine.theta{j,i},2);
    c = copulapdf(vine.families{j,i},[Fp(:,l,j) Fp(:,k,i)],vine.theta{j,i});
    logf(:,i,j) = log(c) + logf(:,l,j);
elseif ~vine.margins{i}.iscont && vine.margins{j}.iscont
    C1 = copulacdf(vine.families{j,i},[Fp(:,l,j) Fp(:,k,i)],vine.theta{j,i});
    C2 = copulacdf(vine.families{j,i},[Fp(:,l,j) Fm(:,k,i)],vine.theta{j,i});
    Fp(:,i,j) = real(exp(log(C1 - C2) - logf(:,k,i)));
    c1 = copulaccdf(vine.families{j,i},[Fp(:,l,j) Fp(:,k,i)],vine.theta{j,i},1);
    c2 = copulaccdf(vine.families{j,i},[Fp(:,l,j) Fm(:,k,i)],vine.theta{j,i},1);
    logf(:,i,j) = log(c1 - c2) + logf(:,l,j) - logf(:,k,i);
elseif vine.margins{i}.iscont && ~vine.margins{j}.iscont
    Fp(:,i,j) = copulaccdf(vine.families{j,i},[Fp(:,l,j) Fp(:,k,i)],vine.theta{j,i},2);
    Fm(:,i,j) = copulaccdf(vine.families{j,i},[Fm(:,l,j) Fp(:,k,i)],vine.theta{j,i},2);
    logf(:,i,j) = log(Fp(:,i,j) - Fm(:,i,j));
else
    C1 = copulacdf(vine.families{j,i},[Fp(:,l,j) Fp(:,k,i)],vine.theta{j,i});
    C2 = copulacdf(vine.families{j,i},[Fp(:,l,j) Fm(:,k,i)],vine.theta{j,i});
    Fp(:,i,j) = real(exp(log(C1 - C2) - logf(:,k,i)));
    C1 = copulacdf(vine.families{j,i},[Fm(:,l,j) Fp(:,k,i)],vine.theta{j,i});
    C2 = copulacdf(vine.families{j,i},[Fm(:,l,j) Fm(:,k,i)],vine.theta{j,i});
    Fm(:,i,j) = real(exp(log(C1 - C2) - logf(:,k,i)));
    logf(:,i,j) = log(Fp(:,i,j) - Fm(:,i,j));
end

end


function evalleftnode(vine,i,j,k,l)

global Fp Fm logf;

if vine.margins{i}.iscont && vine.margins{j}.iscont
    Fp(:,i,j) = copulaccdf(vine.families{i,j},[Fp(:,k,i) Fp(:,l,j)],vine.theta{i,j},1);
    c = copulapdf(vine.families{i,j},[Fp(:,k,i) Fp(:,l,j)],vine.theta{i,j});
    logf(:,i,j) = log(c) + logf(:,l,j);
elseif ~vine.margins{i}.iscont && vine.margins{j}.iscont
    C1 = copulacdf(vine.families{i,j},[Fp(:,k,i) Fp(:,l,j)],vine.theta{i,j});
    C2 = copulacdf(vine.families{i,j},[Fm(:,k,i) Fp(:,l,j)],vine.theta{i,j});
    Fp(:,i,j) = real(exp(log(C1 - C2) - logf(:,k,i)));
    c1 = copulaccdf(vine.families{i,j},[Fp(:,k,i) Fp(:,l,j)],vine.theta{i,j},2);
    c2 = copulaccdf(vine.families{i,j},[Fm(:,k,i) Fp(:,l,j)],vine.theta{i,j},2);
    logf(:,i,j) = log(c1 - c2) + logf(:,l,j) - logf(:,k,i);
elseif vine.margins{i}.iscont && ~vine.margins{j}.iscont
    Fp(:,i,j) = copulaccdf(vine.families{i,j},[Fp(:,k,i) Fp(:,l,j)],vine.theta{i,j},1);
    Fm(:,i,j) = copulaccdf(vine.families{i,j},[Fp(:,k,i) Fm(:,l,j)],vine.theta{i,j},1);
    logf(:,i,j) = log(Fp(:,i,j) - Fm(:,i,j));
else
    C1 = copulacdf(vine.families{i,j},[Fp(:,k,i) Fp(:,l,j)],vine.theta{i,j});
    C2 = copulacdf(vine.families{i,j},[Fm(:,k,i) Fp(:,l,j)],vine.theta{i,j});
    Fp(:,i,j) = real(exp(log(C1 - C2) - logf(:,k,i)));
    C1 = copulacdf(vine.families{i,j},[Fp(:,k,i) Fm(:,l,j)],vine.theta{i,j});
    C2 = copulacdf(vine.families{i,j},[Fm(:,k,i) Fm(:,l,j)],vine.theta{i,j});
    Fm(:,i,j) = real(exp(log(C1 - C2) - logf(:,k,i)));
    logf(:,i,j) = log(Fp(:,i,j) - Fm(:,i,j));
end

end


function evalctree(vine,i,j)

global newnode;

if i == 1
    % Leaf
    evalleftnode(vine,1,j,1,j);
else
    % Node
    if newnode(i-1,i)
        evalctree(vine,i-1,i);
    end
    if newnode(i-1,j)
        evalctree(vine,i-1,j);
    end
    evalleftnode(vine,i,j,i-1,i-1);
end
newnode(i,j) = false;

end
