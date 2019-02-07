function t = dassert(A,B,tol)

% This function tests the equality of two objects.

% Copyright (C) 2011-2017 Dynare Team
%
% This file is part of Dynare (m-unit-tests module).
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare's m-unit-tests module is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
% or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
% more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if ( (nargin<3) || isempty(tol) )
    use_isequal_matlab_builtin = 1;
else
    use_isequal_matlab_builtin = 0;
end

% Shape.
nA = ndims(A);
nB = ndims(B);

% Dimensions.
sA = size(A);
sB = size(B);

% Object type.
cA = class(A);
cB = class(B);

% Number of elements.
nn = prod(sA);

if ~(isequal(nA,nB) && isequal(sA,sB))
    error('dassert:: Dimensions of compared objects A and B don''t match!')
end

if ~isequal(cA,cB)
    error('dassert:: Compared objects are not of the same type!')
end

if isa(A,'double')
    if use_isequal_matlab_builtin
        t = isequal(A,B);
        if ~t
            if ~isoctave && matlab_ver_less_than('7.14')
                t = isequalwithequalnans(A,B);
            else
                t = isequaln(A,B);
            end
        end
    else
        if max(abs(A(:)-B(:)))>tol
            t = 0;
        else
            t = 1;
        end
    end
elseif isa(A,'char')
    t = isequal(A, B);
elseif isa(A,'cell')
    rA = reshape(A, prod(sA), 1);
    rB = reshape(B, prod(sB), 1);
    for i=1:nn
        if nargin>2
            t = dassert(rA{i}, rB{i}, tol);
        else
            t = dassert(rA{i}, rB{i});
        end
        if ~t
            break
        end
    end
elseif isa(A,'struct')
    A = struct2cell(A);
    B = struct2cell(B);
    if nargin>2
        t = dassert(A, B, tol);
    else
        t = dassert(A, B);
    end
else
    if use_isequal_matlab_builtin
        t = isequal(A, B);
    else
        t = isequal(A, B, tol);
    end
end