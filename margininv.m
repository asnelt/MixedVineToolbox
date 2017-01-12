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

function x = margininv(margin,c)
% MARGININV Inverse of univariate margin cumulative distribution function.
%   X = MARGININV(MARGIN,C) returns the inverse of the cumulative
%   distribution function at the values in C of the univariate distribution
%   specified in the struct MARGIN. The margin distribution is specified in
%   the field MARGIN.dist which can be one of
%    'norm'  for the Gaussian distribution,
%    'gam'   for the gamma distribution,
%    'poiss' for the Poisson distribution,
%    'bino'  for the binomial distribution, and
%    'nbin'  for the negative binomial distribution.
%   The parameter of the distribution is specified in the field
%   MARGIN.theta. The size of this field depends on the distribution.
%   C is a vector of length N, where N is the number of samples.
%
%   X is an N-dimensional vector.

% Argument checks
if nargin < 2
    error('margininv: Usage x = margininv(margin,c)');
end
if ~isstruct(margin)
    error('margininv: Argument "margin" must be a struct');
end
if ~isvector(c)
    error('margininv: Argument "c" must be a vector');
end

switch lower(margin.dist)
    case 'norm'
        x = norminv(c,margin.theta(1),margin.theta(2));
    case 'gam'
        x = gaminv(c,margin.theta(1),margin.theta(2));
    case 'poiss'
        x = poissinv(c,margin.theta);
    case 'bino'
        x = binoinv(c,margin.theta(1),margin.theta(2));
    case 'nbin'
        x = nbininv(c,margin.theta(1),margin.theta(2));
    otherwise
        error(['margininv: Unknown margin distribution "' margin.dist '"']);
end

end
