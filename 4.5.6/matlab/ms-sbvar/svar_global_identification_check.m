function indent = svar_global_identification_check(options_)
% function svar_global_identification_check(options_.ms) checks
% identification of s structural VAR
%
% INPUTS
%    options_ms:    (struct)    options
%
% OUTPUTS
%    ident:          (boolean) false = not identified; true = identified
%
% SPECIAL REQUIREMENTS
%    none
%
% REFERENCES
%   J. Rubio Ramirez, D. Waggoner, T. Zha (2010) "Structural Vector
%   Autoregressions: Theory of Identification and Algorithms for
%   Inference" in Review of Economic Studies, 77, 665-696.

% Copyright (C) 2015-2017 Dynare Team
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

ident = false;

if isequal(options_.ms.restriction_fname, 'upper_cholesky') || ...
        isequal(options_.ms.restriction_fname, 'lower_cholesky')
    ident = true;
    if ~options_.noprint
        disp(' ')
        disp('SBVAR: Cholesky identification is always identified')
        disp(' ')
    end
    return
end
nvar = length(options_.varobs);   % number of endogenous variables
nexo = 1;

[Uiconst,Viconst,n0,np,ixmC0Pres,Qi,Ri] = exclusions(nvar,nexo,options_.ms );

% order column constraints by rank
QQranks = zeros(nvar,2);
for j=1:nvar
    n = nvar*(options_.ms.nlags+1);
    QQi{j} = zeros(n,n);
    QQi{j}(1:nvar,1:nvar) = Qi{j};
    QQi{j}(nvar+1:end,nvar+1:end) = Ri{j}(1:end-1,1:end-1);
    QQranks(j,:) = [j,rank(QQi{j})];
end
QQranks = sortrows(QQranks,-2);

ident = true;

for j=1:nvar
    i = QQranks(j,1);
    for k=1:1
        M = [QQi{i}*rand(size(QQi{i},1),nvar);[eye(j) ...
                            zeros(j,nvar-j)]];
        if rank(M) < nvar
            ident = false;
            break
        end
    end
    if ~ident
        break
    end
end

if ~options_.noprint
    disp(' ')
    if ident
        disp('The sufficient condition for SBVAR identification is met')
    else
        disp('WARNGING: The sufficient condition for SBVAR identification is not met')
    end
    disp(' ')
end