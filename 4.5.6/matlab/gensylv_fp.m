function X = gensylv_fp(A, B, C, D, block, tol)
% function X = gensylv_fp(A, B, C, D)
% Solve the Sylvester equation:
% A * X + B * X * C + D = 0
% INPUTS
%   A
%   B
%   C
%   D
%   block : block number (for storage purpose)
%   tol : convergence criteria
% OUTPUTS
%   X solution
%
% ALGORITHM
%   fixed point method
%   MARLLINY MONSALVE (2008): "Block linear method for large scale
%   Sylvester equations", Computational & Applied Mathematics, Vol 27, n°1,
%   p47-59
%   ||A^-1||.||B||.||C|| < 1 is a suffisant condition:
%    - to get a unique solution for the Sylvester equation
%    - to get a convergent fixed-point algorithm
%
% SPECIAL REQUIREMENTS
%   none.
% Copyright (C) 1996-2017 Dynare Team
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


%tol = 1e-12;
evol = 100;
A1 = inv(A);
eval(['persistent hxo_' int2str(block) ';']);
hxo = eval(['hxo_' int2str(block) ';']);
if isempty(hxo)
    X = zeros(size(B, 2), size(C, 1));
else
    X = hxo;
end
it_fp = 0;
maxit_fp = 1000;
Z = - (B * X * C + D);
while it_fp < maxit_fp && evol > tol
    %X_old = X;
    %X = - A1 * ( B * X * C + D);
    %evol = max(max(abs(X - X_old)));
    X = A1 * Z;
    Z_old = Z;
    Z = - (B * X * C + D);
    evol = max(sum(abs(Z - Z_old))); %norm_1
                                     %evol = max(sum(abs(Z - Z_old)')); %norm_inf
    it_fp = it_fp + 1;
end
%fprintf('sylvester it_fp=%d evol=%g | ',it_fp,evol);
if evol < tol
    eval(['hxo_' int2str(block) ' = X;']);
else
    error(['convergence not achieved in fixed point solution of Sylvester equation after ' int2str(it_fp) ' iterations']);
end