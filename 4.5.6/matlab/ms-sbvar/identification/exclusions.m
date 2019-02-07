function [Ui,Vi,n0,np,ixmC0Pres,Qi,Ri] = exclusions(nvar,nexo,options_ms)
% function [Ui,Vi,n0,np,ixmC0Pres] = exclusions(nvar,nexo,options_ms)
%
% INPUTS
%    nvar:                      number endogenous variables
%    nexo:                      number exogenous variables
%    options_ms:    (struct)    options
%
% OUTPUTS
%    Ui
%    Vi
%    n0
%    np
%    ixmC0Pres
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2011-2017 Dynare Team
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

indxC0Pres = options_ms.cross_restrictions;
nlags = options_ms.nlags;

Qi1 = options_ms.Qi;
Ri1 = options_ms.Ri;

Ui = cell(nvar,1);
Vi = cell(nvar,1);
n0 = zeros(nvar,1);
np = zeros(nvar,1);

k = nlags*nvar+1;

for n=1:nvar
    Qi{n} = zeros(nvar,nvar);
    sQ = size(Qi1{n});
    if all(sQ) > 0
        Qi{n}(1:sQ(1),1:sQ(2)) = Qi1{n};
    end
    Ri{n} = zeros(k,k);
    sR = size(Ri1{n});
    if all(sR) > 0
        Ri{n}(1:sR(1),1:sR(2)) = Ri1{n};
    end

    if options_ms.constants_exclusion
        Ri{n}(sR(1)+1,k) = 1;
    end

    Ui{n} = null(Qi{n});
    Vi{n} = null(Ri{n});
    n0(n) = size(Ui{n},2);
    np(n) = size(Vi{n},2);
end

ixmC0Pres = NaN;