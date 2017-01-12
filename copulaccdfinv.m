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

function c = copulaccdfinv(family,u,theta,dim)
% COPULACCDFINV Inverse of copula conditional CDF.
%   C = COPULACCDFINV(FAMILY,U,THETA,DIM) returns the inverse of the copula
%   conditional cumulative distribution function
%   (\partial C) / (\partial U_DIM) of the bivariate copula C with family
%   FAMILY and parameter THETA at the values in U. U has the size [N, 2],
%   where N is the number of samples. FAMILY can be one of
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
    error('copulaccdfinv: Usage c = copulaccdfinv(family,u,theta,dim)');
end
if ~ischar(family)
    error('copulaccdfinv: Argument "family" must be a string');
end
if ~ismatrix(u)
    error('copulaccdfinv: Argument "u" must be a matrix');
end
if size(u,2) ~= 2
    error('copulaccdfinv: Second dimension of u must be 2');
end
if ~ismatrix(theta)
    error('copulaccdfinv: Argument "theta" must be a matrix');
end
if nargin < 4
    dim = 2;
elseif ~isscalar(dim) || (dim ~= 1 && dim ~= 2)
    error('copulaccdfinv: dim must be 1 or 2');
end

family = lower(family);
% Default is dim == 2; if dim == 1, change rotation and arguments
% accordingly
if (dim == 1)
    if strcmp(family,'claytonrot090')
        family = 'claytonrot270';
    elseif strcmp(family,'claytonrot270')
        family = 'claytonrot090';
    end
    u = [u(:,2) u(:,1)];
end

u(u < 0) = 0;
u(u > 1) = 1;

switch family
    case 'ind'
        c = u(:,1);
    case 'gaussian'
        if ~isscalar(theta) || theta < -1 || theta > 1
            error('copulaccdfinv: For Gaussian, theta must be in [-1, 1]');
        end
        x = norminv(u);
        c = normcdf(x(:,1) .* sqrt(1-theta.^2) + theta .* x(:,2));
    case 'student'
        if numel(theta) ~= 2
            error('copulaccdfinv: For Student, theta must be a vector with two elements');
        end
        if theta(1) < -1 || theta(1) > 1
            error('copulaccdfinv: For Student, theta(1) must be in [-1, 1]');
        end
        if theta(2) <= 0
            error('copulaccdfinv: For Student, theta(2) must be positive');
        end
        x = tinv(u,theta(2));
        c = tcdf(sqrt(((1-theta(1).^2) .* (theta(2)+x(:,2).^2)) ./ (theta(2)+1)) .* tinv(u(:,1), theta(2)+1) + theta(1) .* x(:,2),theta(2));
    case 'clayton'
        if ~isscalar(theta) || theta < 0
            error('copulaccdfinv: For Clayton, theta must be in [0, inf)');
        end
        if theta == 0
            c = u(:,1);
        else
            c = (1 - u(:,2).^(-theta) + (u(:,1).*(u(:,2).^(1+theta))).^(-theta./(1+theta))).^(-1./theta);
        end
    case 'claytonrot090'
        if ~isscalar(theta) || theta < 0
            error('copulaccdfinv: For ClaytonRot090, theta must be in [0, inf)');
        end
        if theta == 0
            c = u(:,1);
        else
            c = (1 - (1-u(:,2)).^(-theta) + (u(:,1).*((1-u(:,2)).^(1+theta))).^(-theta./(1+theta))).^(-1./theta);
        end
    case 'claytonrot180'
        if ~isscalar(theta) || theta < 0
            error('copulaccdfinv: For ClaytonRot180, theta must be in [0, inf)');
        end
        if theta == 0
            c = u(:,1);
        else
            c = 1 - (1 - (1-u(:,2)).^(-theta) + ((1-u(:,1)).*((1-u(:,2)).^(1+theta))).^(-theta./(1+theta))).^(-1./theta);
        end
    case 'claytonrot270'
        if ~isscalar(theta) || theta < 0
            error('copulaccdfinv: For ClaytonRot270, theta must be in [0, inf)');
        end
        if theta == 0
            c = u(:,1);
        else
            c = 1 - (1 - u(:,2).^(-theta) + ((1-u(:,1)).*(u(:,2).^(1+theta))).^(-theta./(1+theta))).^(-1./theta);
        end
    otherwise
        error(['copulaccdfinv: Unknown family "' family '"']);
end

% Numerical inaccuracies can cause imaginary parts
c = real(c);

end
