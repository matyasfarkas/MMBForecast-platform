function block = get_internal_doc_block(fname,fpath)
% Extract doc sections from matlab's routine.

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

% Original author: stephane DOT adjemian AT univ DASH lemans DOT fr

% Default output
block = [];

% Open the matlab file.
mid = fopen([fpath '/' fname '.m'],'r');

% Read the matlab file.
file = textscan(mid,'%s','delimiter','\n');
file = file{1};

% Close the matlab file.
fclose(mid);

% Locate the test blocks.
b1 = find(strncmp(file,'%@info:',7))+1;
b2 = find(strncmp(file,'%@eod:',6))-1;
b  = find(strncmp(file,'%!',2));

if ( isempty(b1) && isempty(b2) && isempty(b) )
    % No internal documentation available.
    return
else
    if ( (~isempty(b1) && isempty(b2) && isempty(b)) || ...
         (isempty(b1) && ~isempty(b2) && isempty(b)) || ...
         (isempty(b1) && isempty(b2) && ~isempty(b)) || ...
         (isempty(b1) && ~isempty(b2) && ~isempty(b)) || ...
         (~isempty(b1) && isempty(b2) && ~isempty(b)) || ...
         (~isempty(b1) && ~isempty(b2) && isempty(b)) )
        error('get_internal_doc_block:: There is a problem with the internal block definition!')
    end
    if ( b2~=b(end) || b1~=b(1) || any(b-transpose(b1:1:b2)) )
        error('get_internal_doc_block:: There is a problem with the internal block definition!')
    end
end

% Extract the internal documentation block.
for i=b(1):1:b(end)
    str = file{i};
    block = char(block, str(3:end));
end
block = block(2:end,:);