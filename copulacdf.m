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

function c = copulacdf(family,u,theta)
% COPULACDF Copula cumulative distribution function.
%   C = COPULACDF(FAMILY,U,THETA) returns the copula cumulative
%   distribution function of the bivariate copula family FAMILY with
%   parameter THETA at the values in U. U has the size [N, 2], where N is
%   the number of samples. FAMILY can be one of
%    'ind'           for independence,
%    'gaussian'      for the Gaussian copula family,
%    'student'       for the student copula family,
%    'clayton'       for the Clayton copula family,
%    'claytonrot090' for the 90° clockwise rotated Clayton copula family,
%    'claytonrot180' for survival Clayton copula family, and
%    'claytonrot270' for the 270° clockwise rotated Clayton copula family.
%
%   C is an N-dimensional vector.

% Argument checks
if nargin < 3
    error('copulacdf: Usage c = copulacdf(family,u,theta)');
end
if ~ischar(family)
    error('copulacdf: Argument "family" must be a string');
end
if ~ismatrix(u)
    error('copulacdf: Argument "u" must be a matrix');
end
if size(u,2) ~= 2
    error('copulacdf: Second dimension of u must be 2');
end
if ~ismatrix(theta)
    error('copulacdf: Argument "theta" must be a matrix');
end

family = lower(family);
u(u < 0) = 0;
u(u > 1) = 1;

switch family
    case 'ind'
        c = prod(u,2);
    case 'gaussian'
        if ~isscalar(theta) || theta < -1 || theta > 1
            error('copulacdf: For Gaussian, theta must be in [-1, 1]');
        end
        c = mvncdf(norminv(u),[0 0],[1 theta; theta 1]);
    case 'student'
        if numel(theta) ~= 2
            error('copulacdf: For Student, theta must be a vector with two elements');
        end
        if theta(1) < -1 || theta(1) > 1
            error('copulacdf: For Student, theta(1) must be in [-1, 1]');
        end
        if theta(2) <= 0
            error('copulacdf: For Student, theta(2) must be positive');
        end
        c = mvtcdf(tinv(u,theta(2)),[1 theta(1); theta(1) 1],theta(2));
    case 'clayton'
        if ~isscalar(theta) || theta < 0
            error('copulacdf: For Clayton, theta must be in [0, inf)');
        end
        if theta == 0
            c = prod(u,2);
        else
            c = max(u(:,1).^(-theta) + u(:,2).^(-theta) - 1,0).^(-1./theta);
        end
    case 'claytonrot090'
        if ~isscalar(theta) || theta < 0
            error('copulacdf: For ClaytonRot090, theta must be in [0, inf)');
        end
        if theta == 0
            c = prod(u,2);
        else
            c = u(:,1) - max(u(:,1).^(-theta) + (1-u(:,2)).^(-theta) - 1,0).^(-1./theta);
        end
    case 'claytonrot180'
        if ~isscalar(theta) || theta < 0
            error('copulacdf: For ClaytonRot180, theta must be in [0, inf)');
        end
        if theta == 0
            c = prod(u,2);
        else
            c = u(:,1) + u(:,2) - 1 + max((1-u(:,1)).^(-theta) + (1-u(:,2)).^(-theta) - 1,0).^(-1./theta);
        end
    case 'claytonrot270'
        if ~isscalar(theta) || theta < 0
            error('copulacdf: For ClaytonRot270, theta must be in [0, inf)');
        end
        if theta == 0
            c = prod(u,2);
        else
            c = u(:,2) - max((1-u(:,1)).^(-theta) + u(:,2).^(-theta) - 1,0).^(-1./theta);
        end
    otherwise
        error(['copulacdf: Unknown family "' family '"']);
end

end
