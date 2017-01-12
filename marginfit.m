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

function margin = marginfit(x,iscont)
% MARGINFIT Univariate margin estimates.
%   MARGIN = MARGINFIT(X,ISCONT) selects the best fitting distribution
%   MARGIN.dist given the data U, where U has the size [N, 2] for N the
%   number of samples. The selection criterion is the Akaike information
%   criterion. The selected margin is specified in the field MARGIN.dist as
%   one of
%    'norm'  for the Gaussian distribution or
%    'gam'   for the gamma distribution
%   if ISCONT is true, or one of
%    'poiss' for the Poisson distribution,
%    'bino'  for the binomial distribution, and
%    'nbin'  for the negative binomial distribution
%   if ISCONT is false. The maximum likelihood estimate (MLE) of the
%   parameter MARGIN.theta of the distribution given the data U is
%   specified in the field MARGIN.theta.
%
%   The size of MARGIN.theta depends on the selected distribution
%   MARGIN.dist. MARGIN.iscont is set to ISCONT.

% Argument checks
if nargin < 2
    error('marginfit: Usage margin = marginfit(x,iscont)');
end
if ~isvector(x)
    error('marginfit: Argument "x" must be a vector');
end
if ~islogical(iscont) || ~isscalar(iscont)
    error('marginfit: Argument "iscont" must be a boolean');
end

if iscont
    % Continuous margin
    if any(x <= 0)
        dists = {'norm'};
    else
        dists = {'norm','gam'};
    end
else
    % Discrete margin
    if mean(x) >= var(x)
        dists = {'poiss','bino'};
    else
        dists = {'poiss','bino','nbin'};
    end
end

aic = zeros(length(dists),1);
theta = cell(length(dists),1);
for i = 1:length(dists)
    [logp,theta{i}] = logpdf(dists{i},x);
    aic(i) = 2*length(theta{i}) - 2*logp;
end
% Use AIC for selecting the best distribution
[~,imin] = min(aic);
margin.dist = dists{imin};
margin.theta = theta{imin};
margin.iscont = iscont;

end


function [logp,theta] = logpdf(dist,x)

switch dist
    case 'norm'
        [theta(1),theta(2)] = normfit(x);
    case 'gam'
        theta = gamfit(x);
    case 'poiss'
        theta = poissfit(x);
    case 'bino'
        theta(1) = max(x);
        theta(2) = mean(binofit(x,theta(1)));
    case 'nbin'
        theta = nbinfit(x);
    otherwise
        error(['marginfit: Unknown distribution "' dist '"']);
end
margin.dist = dist;
margin.theta = theta;
logp = sum(log(marginpdf(margin,x)));

end
