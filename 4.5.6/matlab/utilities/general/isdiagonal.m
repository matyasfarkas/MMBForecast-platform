function b = isdiagonal(A) % --*-- Unitary tests --*--

% Copyright (C) 2014-2017 Dynare Team
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

if isnumeric(A)
    if isquare(A)
        % Find non zero elements in matrix A...
        [ir, ic] = find(A);
        % If the non zero elements are on the diagonal, the corresponding elements
        % in ir and ic (row and column numbers) should be equal.
        b = isequal(ir, ic);
    else
        error('isdiagonal: Input must be a square matrix!')
    end
else
    error('isdiagonal: Input must be numeric!')
end

%@test:1
%$ A = zeros(3,3);
%$ t = isdiagonal(A);
%$ T = all(t);
%@eof:1

%@test:2
%$ A = zeros(3,3); A(1,3) = 1;
%$ t = ~isdiagonal(A);
%$ T = all(t);
%@eof:2

%@test:3
%$ A = randn(3,3);
%$ t = ~isdiagonal(A);
%$ T = all(t);
%@eof:3

%@test:4
%$ A = diag(randn(3,1));
%$ t = isdiagonal(A);
%$ T = all(t);
%@eof:4

%@test:5
%$ A = diag(randn(3,1)); A(1,1) = 0;
%$ t = isdiagonal(A);
%$ T = all(t);
%@eof:5
