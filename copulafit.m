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

function theta = copulafit(family,u)
% COPULAFIT Copula parameter estimates.
%   THETA = COPULAFIT(FAMILY,U) returns the maximum likelihood estimate
%   (MLE) of the parameter THETA of the bivariate copula family FAMILY,
%   given the data U, where U has the size [N, 2] for N the number of
%   samples. FAMILY can be one of
%    'ind'           for independence,
%    'gaussian'      for the Gaussian copula family,
%    'student'       for the student copula family,
%    'clayton'       for the Clayton copula family,
%    'claytonrot090' for the 90° clockwise rotated Clayton copula family,
%    'claytonrot180' for survival Clayton copula family, and
%    'claytonrot270' for the 270° clockwise rotated Clayton copula family.
%
%   The size of THETA depends on the copula family FAMILY.

% Argument checks
if nargin < 2
    error('copulafit: Usage theta = copulafit(family,u)');
end
if ~ischar(family)
    error('copulafit: Argument "family" must be a string');
end
if ~ismatrix(u)
    error('copulafit: Argument "u" must be a matrix');
end
if size(u,2) ~= 2
    error('copulafit: Second dimension of u must be 2');
end

family = lower(family);
u(u < 0) = 0;
u(u > 1) = 1;

switch family
    case 'ind'
        theta = [];
    case 'gaussian'
        fun = @(theta) -sum(log(copulapdf(family,u,theta)+eps));
        theta = fminbnd(fun,-1,1,optimset('Display','off'));
    case 'student'
        fun = @(theta) -sum(log(copulapdf(family,u,theta)+eps));
        init = [0; 1];
        lb = [-1; 1e-1];
        ub = [1; 1000];
        options = optimset('Algorithm','interior-point','Display','off');
        theta = fmincon(fun,init,[],[],[],[],lb,ub,[],options);
    case {'clayton','claytonrot090','claytonrot180','claytonrot270'}
        fun = @(theta) -sum(log(copulapdf(family,u,theta)+eps));
        theta = fminbnd(fun,1e-3,20,optimset('Display','off'));
    otherwise
        error(['copulafit: Unknown family "' family '"']);
end

end
