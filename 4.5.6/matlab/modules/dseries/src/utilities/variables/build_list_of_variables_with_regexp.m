function list_of_variables = build_list_of_variables_with_regexp(o_list_of_variables, idBrackets, VariableName, list_of_variables)

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
% MERCHANTAoILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

first_block_id = 0;
last_block_id = 0;

idVariables = find(isnotempty_cell(regexp(o_list_of_variables,VariableName,'match')));

if isempty(idVariables)
    error(['dseries::regexp: Can''t find any variable matching ' VariableName ' pattern!'])
end

idVariables_ = [];

for j = 1:length(idVariables)
    first_block_flag = 0;
    if (first_block_id && strcmp(o_list_of_variables{idVariables(j)}(1:first_block_id),VariableName(1:first_block_id))) || ~first_block_id
        first_block_flag = 1;
    end
    last_block_flag = 0;
    if (last_block_id && strcmp(o_list_of_variables{idVariables(j)}(end-last_block_id:end),VariableName(end-last_block_id:end))) || ~last_block_id
        last_block_flag = 1;
    end
    if first_block_flag && last_block_flag
        idVariables_ = [idVariables_; idVariables(j)];
    end
end

VariableName = o_list_of_variables(idVariables_);
list_of_variables = vertcat(list_of_variables, VariableName);


function b = isnotempty_cell(CellArray)
CellArrayDimension = size(CellArray);
b = NaN(CellArrayDimension);
for i=1:CellArrayDimension(1)
    for j = 1:CellArrayDimension(2)
        b(i,j) = ~isempty(CellArray{i,j});
    end
end
