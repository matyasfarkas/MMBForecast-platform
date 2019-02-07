function k = symmetric_matrix_index(i,j,n)
% function k = symmetric_matrix_index(i,j,n)
% Returns index number of variable combination (i,j) in vech(A) where A is
% an symmetric n by n matrix and vech creates row vector by stacking rows
% of A on and above the diagonal
%
% Inputs:
%   i   [scalar]    index of first variable
%   j   [scalar]    index of second variable
%   n   [scalar]    number of variables
% Outputs:
%   k   [scalar]    index of variable combination in vech(A)

% Copyright (C) 2007-2017 Dynare Team
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

k = (i-1)*n+j-i*(i-1)/2;