function o = write(o, fid)
%function o = write(o, fid)
% Write a Vspace object
%
% INPUTS
%   o           [vspace]   vspace object
%   fid         [integer]  file id
%
% OUTPUTS
%   o           [vspace]   vspace object
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2013-2015 Dynare Team
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

assert(fid ~= -1);

for i=1:o.number
    fprintf(fid, ' \\par \\medskip ');
end

if o.hline > 0
    fprintf(fid, '\\\\\n');
    for i=1:o.hline
        fprintf(fid, '\\midrule');
    end
end
end