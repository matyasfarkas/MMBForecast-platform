function dxp=getPowerDeriv(x,p,k)
%function dxp=getPowerDeriv(x,p,k)
% The k-th derivative of x^p
%
% INPUTS
%    x: base
%    p: power
%    k: derivative order
%
% OUTPUTS
%    dxp: k-th derivative of x^p
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2011-2012 Dynare Team
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

if (abs(x) < 1e-12) && (p > 0) && (k > p) && (abs(p - round(p)) < 1e-12)
    dxp = 0;
else
    dxp = x^(p-k);
    for i=0:k-1
        dxp = dxp*p;
        p = p-1;
    end
end
end
