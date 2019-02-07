function [dr,info]=AIM_first_order_solver(jacobia,M,dr,qz_criterium)

%@info:
%! @deftypefn {Function File} {[@var{dr},@var{info}] =} AIM_first_order_solver (@var{jacobia},@var{M},@var{dr},@var{qz_criterium})
%! @anchor{AIM_first_order_solver}
%! @sp 1
%! Computes the first order reduced form of the DSGE model using AIM.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item jacobia
%! Matrix containing the Jacobian of the model
%! @item M
%! Matlab's structure describing the model (initialized by @code{dynare}).
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @item qz_criterium
%! Double containing the criterium to separate explosive from stable eigenvalues
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @item info
%! Integer scalar, error code.
%! @sp 1
%! @table @ @code
%! @item info==0
%! No error.
%! @item info==102
%! roots not correctly computed by real_schur
%! @item info==103
%! Blanchard & Kahn conditions are not satisfied: no stable equilibrium.
%! @item info==104
%! Blanchard & Kahn conditions are not satisfied: indeterminacy.
%! @item info==135
%! too many explosive roots and q(:,right) is singular
%! @item info==145
%! too few big roots, and q(:,right) is singular
%! @item info==105
%! q(:,right) is singular
%! @item info==161
%! too many exact siftrights
%! @item info==162
%! too many numeric shiftrights
%! @end table
%! @end table
%! @end deftypefn
%@eod:

% Copyright (C) 2001-2017 Dynare Team
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

info = 0;

[dr,aimcode]=dynAIMsolver1(jacobia,M,dr);

if aimcode ~=1
    info(1) = convertAimCodeToInfo(aimCode); %convert to be in the 100 range
    info(2) = 1.0e+8;
    return
end
A = kalman_transition_matrix(dr,M.nstatic+(1:M.nspred), 1:M.nspred,...
                             M.exo_nbr);
dr.eigval = eig(A);
disp(dr.eigval)
nd = size(dr.kstate,1);
nba = nd-sum( abs(dr.eigval) < qz_criterium );

nsfwrd = M.nsfwrd;

if nba ~= nsfwrd
    temp = sort(abs(dr.eigval));
    if nba > nsfwrd
        temp = temp(nd-nba+1:nd-nsfwrd)-1-qz_criterium;
        info(1) = 3;
    elseif nba < nsfwrd
        temp = temp(nd-nsfwrd+1:nd-nba)-1-qz_criterium;
        info(1) = 4;
    end
    info(2) = temp'*temp;
    return
end
