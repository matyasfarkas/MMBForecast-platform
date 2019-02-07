function c = read_key_value_string(s)

% Transforms a string of Key-Value options as returned by the preprocessor (for optim option in the
% estimation command) into a cell (first column for the option name ans second column for the
% option value).

% Copyright (C) 2011-2017 Dynare Team
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

i_opening_bracket = strfind(s,'(');
i_closing_bracket = strfind(s,')');

if ~isequal(length(i_opening_bracket),length(i_closing_bracket))
    error('read_key_value_string: The number of opening and closing brackets does not match!')
end

%deal with sublists
i_begin_sublist = strfind(s,'''(');
i_end_sublist = strfind(s,')''');

if ~isequal(length(i_begin_sublist),length(i_end_sublist))
    error('read_key_value_string: Sublists could not be identified!')
end

iComma = strfind(s,',');

%delete commata in sublists from further checks
for sublist_iter=length(i_begin_sublist):-1:1
    iComma(iComma>=i_begin_sublist(sublist_iter) & iComma<=i_end_sublist(sublist_iter))=[];
end

nComma = length(iComma);

if iseven(nComma)
    error('read_key_value_string: Wrong number of Key-Value pairs!')
end

c = cell((nComma+1)/2,2);

for i = 1:nComma
    j = comma2opt(i);
    if j>0
        continue
    end
    if isequal(i,1)
        i1 = 1;
        i2 = iComma(i)-1;
    else
        i1 = iComma(i-1)+1;
        i2 = iComma(i)-1;
    end
    if isequal(i,nComma)
        i3 = iComma(i)+1;
        i4 = length(s);
    else
        i3 = iComma(i)+1;
        i4 = iComma(i+1)-1;
    end
    c(-j,1) = {s(i1+1:i2-1)};
    tmp = str2num(s(i3:i4));
    if isempty(tmp)
        if strcmp(s(i3:i3+1),'''(') && strcmp(s(i4-1:i4),')''') %sublist, delete brackets
            c(-j,2) = {s(i3+2:i4-2)};
        else
            c(-j,2) = {s(i3+1:i4-1)};
        end
    else
        c(-j,2) = {tmp};
    end
end


function j = comma2opt(i)
if isodd(i)
    % The comma is a separator between a Key and a Value (returned j is minus the option number).
    j = - (i+1)/2;
else
    % The comma is a separator between two options (returned j is the option number).
    j = i/2;
end