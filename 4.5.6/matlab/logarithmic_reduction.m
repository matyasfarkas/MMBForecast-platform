function [X1, info] = logarithmic_reduction(A,B,C,tol,maxit,check)

%@info:
%! @deftypefn {Function File} {[@var{X1}, @var{info}] =} logarithmic_reduction (@var{A},@var{B},@var{C},@var{tol},@var{maxit},@var{check})
%! @anchor{logarithmic_reduction}
%! @sp 1
%! Solves the quadratic matrix equation AX^2 + BX + C = 0.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A
%! Square matrix of doubles, n*n.
%! @item B
%! Square matrix of doubles, n*n.
%! @item C
%! Square matrix of doubles, n*n.
%! @item tol
%! Scalar double, tolerance parameter.
%! @item maxit
%! Scalar integer, maximum number of iterations.
%! @item check
%! Scalar integer, if non zero the solution is checked.
%! @end table
%! @sp 1
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item X1
%! Square matrix of doubles, n*n, solution of the matrix equation.
%! @item info
%! Scalar integer, if nonzero the algorithm failed in finding the solution of the matrix equation.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
%! @strong{References:}
%! @sp 1
%! G. Latouche and V. Ramaswami (1993), "A logarithmic reduction algorithm for Quasi-Birth-Death processes", in Journal of Applied Probability, Vol. 30, No. 3, pp. 650-674.
%! @sp 2
%! @end deftypefn
%@eod:

% Copyright (C) 2012 Dynare Team
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

n = length(A);
i = 1:n;

tmp0 = -B\[A,C];

X0 = tmp0(:,n+i);
U0 = tmp0(:,i);

kk = 0;
cc = 1;

while cc>tol && kk<=maxit
    tmp1 = (eye(n)-tmp0*[tmp0(:,n+i);tmp0(:,i)])\[tmp0(:,i)*tmp0(:,i),tmp0(:,n+i)*tmp0(:,n+i)];
    X1 = X0 + U0*tmp1(:,n+i);
    U1 = U0*tmp1(:,i);
    cc = max(max(abs(X1-X0)));
    X0 = X1; U0 = U1;
    tmp0 = tmp1;
    kk = kk+1;
end

if kk==maxit
    disp(['logarithmic_reduction:: Convergence not achieved after ' int2str(maxit) ' iterations!']);
    info = 1;
end

if nargin>5 && check
    if max(max(abs(A*X1*X1 + B*X1 + C)))>tol
        disp(['logarithmic_reduction:: Algotithm did not converge to the solution of the matrix quadratic equation!']);
        info = 1;
    end
end