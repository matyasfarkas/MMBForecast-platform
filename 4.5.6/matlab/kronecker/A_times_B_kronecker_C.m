function [D, err] = A_times_B_kronecker_C(A,B,C,fake)

%@info:
%! @deftypefn {Function File} {[@var{D}, @var{err}] =} A_times_B_kronecker_C (@var{A},@var{B},@var{C},@var{fake})
%! @anchor{kronecker/A_times_B_kronecker_C}
%! @sp 1
%! Computes A*kron(B,C).
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item A
%! mA*nA matrix of doubles.
%! @item B
%! mB*nB matrix of doubles.
%! @item C
%! mC*nC matrix of doubles.
%! @item fake
%! Anything you want, just a fake parameter (because the mex version admits a last argument specifying the number of threads to be used in parallel mode).
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item D
%! mA*(nC*nB) or mA*(nB*nB) matrix of doubles.
%! @item err
%! Integer scalar equal to zero (if all goes well).
%! @end table
%! @sp 2
%! @strong{Remarks}
%! @sp 1
%! [1] This routine is called by Dynare if and only the mex version is not compiled (also used for testing purposes).
%! @sp 1
%! [2] This routine can be called with three or four arguments. In the first case A*kron(B,B) is computed.
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{kronecker/sparse_hessian_times_B_kronecker_C}, @ref{dr1}, @ref{simult_}
%! @sp 2
%! @strong{This function calls:}
%!
%! @end deftypefn
%@eod:

% Copyright (C) 1996-2012 Dynare Team
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

% Original author: stephane DOT adjemian AT univ DASH lemans DOT fr

% Chek number of inputs and outputs.
if nargin>4 || nargin<3
    error('A_times_B_kronecker_C takes 3 or 4 input arguments and provides 2 output arguments.')
end


% Get & check dimensions. Initialization of the output matrix.
[mA,nA] = size(A);
[mB,nB] = size(B);
if nargin == 4
    [mC,nC] = size(C);
    if mB*mC ~= nA
        error('Input dimension error!')
    end
    D = zeros(mA,nB*nC);
    loop = (mB*nB*mC*nC > 1e7);
else
    if mB*mB ~= nA
        error('Input dimension error!')
    end
    D = zeros(mA,nB*nB);
    loop = (mB*nB*mB*nB > 1e7);
end
% Computational part.
if loop
    if nargin == 4
        k1 = 1;
        for i1=1:nB
            for i2=1:nC
                D(:,k1) = A * kron(B(:,i1),C(:,i2));
                k1 = k1 + 1;
            end
        end
    else
        k1 = 1;
        for i1=1:nB
            for i2=1:nB
                D(:,k1) = A * kron(B(:,i1),B(:,i2));
                k1 = k1 + 1;
            end
        end
    end
else
    if nargin == 4
        D = A * kron(B,C);
    else
        D = A * kron(B,B);
    end
end
err = 0;