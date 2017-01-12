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

function x = mixedvinernd(vine,cases)
% MIXEDVINERND Mixed copula vine random numbers.
%   X = MIXEDVINERND(VINE,CASES) generates random numbers from the mixed
%   copula vine specified in the struct VINE. The vine type is specified in
%   the field VINE.type which can only be 'c-vine' for the canonical vine
%   for now.
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
%   X has the size [CASES, D], where CASES is the number of samples and D
%   is the dimension of each sample.

% Argument checks
if nargin < 1
    error('mixedvinernd: Usage x = mixedvinernd(vine,cases)');
end
if ~isstruct(vine)
    error('mixedvinernd: Argument "vine" must be a struct');
end
if nargin < 2
    cases = 1;
elseif ~isscalar(cases)
    error('mixedvinernd: Argument "cases" must be scalar');
end

d = length(vine.margins);
w = rand(cases,d);
v = zeros(cases,d,d);
v(:,1,1) = reshape(w(:,1),[cases 1 1]);
for i = 2:d
    v(:,i,1) = reshape(w(:,i),[cases 1 1]);
    for k = (i-1):-1:1
        v(:,i,1) = copulaccdfinv(vine.families{k,i},[v(:,k,k) v(:,i,1)],vine.theta{k,i},1);
    end
    if i < d
        for j = 1:(i-1)
            v(:,i,j+1) = copulaccdf(vine.families{j,i},[v(:,j,j) v(:,i,j)],vine.theta{j,i},1);
        end
    end
end
u = reshape(v(:,:,1),size(w));
x = zeros(size(u));
for i = 1:d
    x(:,i) = margininv(vine.margins{i},u(:,i));
end

end
