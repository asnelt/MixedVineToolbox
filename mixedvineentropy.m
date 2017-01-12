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

function [h,stderr] = mixedvineentropy(vine,alpha,erreps,cases)
% MIXEDVINEENTROPY Mixed copula vine entropy estimate.
%   [H,STDERR] = MIXEDVINEENTROPY(VINE,ALPHA,ERREPS,CASES) returns an
%   estimate H of the entropy of the mixed copula vine distribution
%   specified in the struct VINE. The vine type is specified in the field
%   VINE.type which can only be 'c-vine' for the canonical vine for now.
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
%   ALPHA is the significance level of the estimate (default ALPHA = 0.05).
%   ERREPS is the maximum standard error of the entropy estimate as a
%   stopping criterion (default ERREPS = 1e-3).
%   CASES is the number of samples that are drawn in each iteration of the
%   Monte Carlo estimation (default CASES = 1000).
%
%   H and STDERR are scalars in unit bit (base 2 logarithm).

% Argument checks
if nargin < 1
    error('mixedvineentropy: Usage [h,stderr] = mixedvineentropy(vine,alpha,erreps,cases)');
end
if ~isstruct(vine)
    error('mixedvineentropy: Argument "vine" must be a struct');
end
if nargin < 2 || isempty(alpha)
    alpha = 0.05;
elseif ~isscalar(alpha)
    error('mixedvineentropy: Argument "alpha" must be a scalar');
end
if nargin < 3 || isempty(erreps)
    erreps = 1e-3;
elseif ~isscalar(erreps)
    error('mixedvineentropy: Argument "erreps" must be a scalar');
end
if nargin < 4
    cases = 1000;
elseif ~isscalar(cases)
    error('mixedvineentropy: Argument "cases" must be a scalar');
end

% Gaussian confidence interval for erreps and level alpha
conf = norminv(1 - alpha,0,1);

stderr = inf;
h = 0;
varsum = 0;
k = 0;
while stderr >= erreps
    % Generate samples
    x = mixedvinernd(vine,cases);
    [~,logp] = mixedvinepdf(vine,x);
    log2p = logp(~isinf(logp)) / log(2);
    k = k + 1;
    % Monte-Carlo estimate of entropy
    h = h + (-mean(log2p) - h) / k;
    % Estimate standard error
    varsum = varsum + sum((-log2p - h) .^ 2);
    stderr = conf * sqrt(varsum / (k * cases * (k * cases - 1)));
end

end
