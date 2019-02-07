function print_moments_implied_prior(ModelInfo, mm, vm, mv, vv)

% This routine prints in the command window some descriptive statistics
% about the endogenous variables implied prior moments.

% Copyright (C) 2016-2017 Dynare Team
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

% First order moments.

disp('(Implied) Prior moments of the endogenous variables'' expectation')
disp(printline(64, '-'))

T1 = 'VARIABLE ';
T2 = sprintf('Prior mean \t Prior st. dev.');

for i=1:ModelInfo.orig_endo_nbr
    Name = deblank(ModelInfo.endo_names(i, :));
    T1 = strvcat(T1, Name);
    str = sprintf(' %6.4f \t %6.4f', mm(i), sqrt(vm(i)));
    T2 = strvcat(T2, str);
end

T0 = repmat('  ', ModelInfo.orig_endo_nbr+1, 1);

TT = [T1, T0, T2];
l0 = printline(size(TT, 2)+1, '-');
TT = strvcat(l0, TT(1,:), l0, TT(2:end,:), l0);

skipline(2)
disp(TT)
skipline(2)

disp('(Implied) Prior moments of the endogenous variables'' variance')
disp(printline(61, '-'))

T1a = 'VARIABLE-1';
T1b = 'VARIABLE-2';
T2a = 'Prior mean';
T2b = 'Prior st.dev.';

for i=1:ModelInfo.orig_endo_nbr
    for j=i:ModelInfo.orig_endo_nbr
        Name1 = deblank(ModelInfo.endo_names(i, :));
        Name2 = deblank(ModelInfo.endo_names(j, :));
        T1a = strvcat(T1a, Name1);
        T1b = strvcat(T1b, Name2);
        sta = sprintf('%12.8f', mv(i,j));
        stb = sprintf('%12.8f', vv(i,j));
        T2a = strvcat(T2a, sta);
        T2b = strvcat(T2b, stb);
    end
end

T0 = repmat('  ', ModelInfo.orig_endo_nbr*(ModelInfo.orig_endo_nbr+1)/2+1, 1);

TT = [T1a, T0, T1b, T0, T2a, T0, T2b];
l0 = printline(size(TT, 2)+1, '-');
TT = strvcat(l0, TT(1,:), l0, TT(2:end,:), l0);

skipline(2)
disp(TT)
skipline(2)