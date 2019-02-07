function  [ldens,Dldens,D2ldens] = lpdfgam(x,a,b)
% Evaluates the logged GAMMA PDF at x.
%
% INPUTS
%    x     [double]  m*n matrix of locations,
%    a     [double]  m*n matrix or scalar, First GAMMA distribution parameters (shape),
%    b     [double]  m*n matrix or scalar, Second GAMMA distribution parameters (scale),
%
% OUTPUTS
%    ldens [double]  m*n matrix of logged GAMMA densities evaluated at x.
%
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2003-2017 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

ldens = -Inf( size(x) ) ;
idx = find( x>0 );

if length(a)==1
    ldens(idx) = -gammaln(a) - a*log(b) + (a-1)*log(x(idx)) - x(idx)/b ;
else
    ldens(idx) = -gammaln(a(idx)) - a(idx).*log(b(idx)) + (a(idx)-1).*log(x(idx)) - x(idx)./b(idx) ;
end



if nargout >1
    if length(a)==1
        Dldens(idx) = (a-1)./(x(idx)) - ones(length(idx),1)/b ;
    else
        Dldens(idx) = (a(idx)-1)./(x(idx)) - ones(length(idx),1)./b(idx) ;
    end
end

if nargout == 3
    if length(a)==1
        D2ldens(idx) = -(a-1)./(x(idx)).^2;
    else
        D2ldens(idx) = -(a(idx)-1)./(x(idx)).^2;
    end
end