function range_indx = get_subsamples_range_indx(subsamples_indx, range_label)
% function range_indx = get_subsamples_range_indx(subsamples_indx, range_label)
%
% Returns the existing estimation_info.subsamples.range_index index
% for the range_label
%
% INPUTS
%    subsamples_indx    [integer] index in estimation_info.subsamples
%    range_label        [string]  label for range
%
% OUTPUTS
%    range_indx         [integer] existing index in the
%                                 estimation_info.subsamples.range_index
%                                 structure associated with the range_label
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2012-2017 Dynare Team
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

global estimation_info

range_indx = find(strcmp(range_label, estimation_info.subsamples(subsamples_indx).range_index) == 1);

if size(range_indx,2) ~= 1
    error(['Error: Index not found in estimation_info.subsamples(' ...
           num2str(subsamples_indx) ').range_index for label ' range_label]);
end
