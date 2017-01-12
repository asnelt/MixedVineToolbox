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

function p = marginpdf(margin,x)
% MARGINPDF Univariate margin probability density function.
%   P = MARGINPDF(MARGIN,X) returns the probability density at the values
%   in X of the univariate distribution specified in the struct MARGIN.
%   The margin distribution is specified in the field MARGIN.dist which can
%   be one of
%    'norm'  for the Gaussian distribution,
%    'gam'   for the gamma distribution,
%    'poiss' for the Poisson distribution,
%    'bino'  for the binomial distribution, and
%    'nbin'  for the negative binomial distribution.
%   The parameter of the distribution is specified in the field
%   MARGIN.theta. The size of this field depends on the distribution.
%   X is a vector of length N, where N is the number of samples.
%
%   P is an N-dimensional vector.

% Argument checks
if nargin < 2
    error('marginpdf: Usage p = marginpdf(margin,x)');
end
if ~isstruct(margin)
    error('marginpdf: Argument "margin" must be a struct');
end
if ~isvector(x)
    error('marginpdf: Argument "x" must be a vector');
end

switch lower(margin.dist)
    case 'norm'
        p = normpdf(x,margin.theta(1),margin.theta(2));
    case 'gam'
        p = gampdf(x,margin.theta(1),margin.theta(2));
    case 'poiss'
        p = poisspdf(x,margin.theta);
    case 'bino'
        p = binopdf(x,margin.theta(1),margin.theta(2));
    case 'nbin'
        p = nbinpdf(x,margin.theta(1),margin.theta(2));
    otherwise
        error(['marginpdf: Unknown margin distribution "' margin.dist '"']);
end

end
