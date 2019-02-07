function [Ui,Vi,n0,np,ixmC0Pres] = upper_cholesky(nvar,nexo,options_ms)
% function [Ui,Vi,n0,np,ixmC0Pres] = upper_cholesky(nvar,nexo,options_ms)
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

lags = options_ms.nlags;
indxC0Pres = options_ms.cross_restrictions;

Ui = cell(nvar,1);
Vi = cell(nvar,1);
n0 = zeros(nvar,1);
np = zeros(nvar,1);

if (nargin==2)
    nexo = 1;
elseif (nargin==3)
    indxC0Pres = 0;
end

k = lags*nvar+nexo;
Qi = zeros(nvar,nvar,nvar);
Ri = zeros(k,k,nvar);

for ii=2:nvar
    Qi(ii-1,ii-1,ii)=1;
    Qi(:,:,ii)=Qi(:,:,ii)+Qi(:,:,ii-1);
end

if options_ms.constants_exclusion
    for i=1:nvar
        Ri(i,k,i) = 1;
    end
end

for n=1:nvar
    Ui{n} = null(Qi(:,:,n));
    Vi{n} = null(Ri(:,:,n));
    n0(n) = size(Ui{n},2);
    np(n) = size(Vi{n},2);
end

ixmC0Pres = NaN;