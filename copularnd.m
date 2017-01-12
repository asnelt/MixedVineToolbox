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

function u = copularnd(family,theta,cases)
% COPULARND Copula random numbers.
%   U = COPULARND(FAMILY,THETA,CASES) generates random numbers from the
%   bivariate copula family FAMILY with parameter THETA. FAMILY can be one
%   of
%    'ind'           for independence,
%    'gaussian'      for the Gaussian copula family,
%    'student'       for the student copula family,
%    'clayton'       for the Clayton copula family,
%    'claytonrot090' for the 90° clockwise rotated Clayton copula family,
%    'claytonrot180' for survival Clayton copula family, and
%    'claytonrot270' for the 270° clockwise rotated Clayton copula family.
%
%   U has the size [CASES, 2], where CASES is the number of bivariate
%   samples.

% Argument checks
if nargin < 2
    error('copularnd: Usage u = copularnd(family,theta,cases)');
end
if ~ischar(family)
    error('copularnd: Argument "family" must be a string');
end
if ~ismatrix(theta)
    error('copularnd: Argument "theta" must be a matrix');
end
if nargin < 3
    cases = 1;
elseif ~isscalar(cases)
    error('copularnd: Argument "cases" must be scalar');
end

family = lower(family);
u = rand(cases,2);
u(:,1) = copulaccdfinv(family,u,theta);

end
