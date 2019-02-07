function [Q, A0, Aplus, Zeta] = load_flat_file(file_tag)
%function [Q, A0, Aplus, Zeta] = load_flat_file(file_tag)
% Loads file saved by save_draws option of ms_simulation.
%
% INPUTS
%   file_tag:    the file tag used with ms_simulation
%
% OUTPUTS
%   Q:           Q matrix
%   A0:          A0 matrix
%   Aplus:       A+ matrix
%   Zeta:        Zeta matrx
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2013 Dynare Team
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

Q = [];
A0 = [];
Aplus = [];
Zeta = [];

try
    headerfid = fopen(['est_flat_header_' file_tag '.out'], 'r');
catch
    error(['Can''t find est_flat_header_' file_tag '.out'])
end
headerfile=textscan(headerfid, '%s');
fclose(headerfid);
flatfile = dlmread(['est_flat_' file_tag '.out'], ' ');

headerfile = headerfile{:};
for i=1:length(headerfile)
    line = char(headerfile{i});
    indob = strfind(line, '[');
    indcb = strfind(line, ']');
    indop = strfind(line, '(');
    indcp = strfind(line, ')');
    indc  = strfind(line, ',');
    if isempty(indob)
        name = line(1:indop-1);
        dim3 = [];
    else
        name = line(1:indob-1);
        dim3 = line(indob+1:indcb-1);
    end
    row = line(indop+1:indc-1);
    col = line(indc+1:indcp-1);

    if isempty(dim3)
        eval([name '(' row ',' col ') = ' num2str(flatfile(i)) ';']);
    else
        eval([name '(' row ',' col ',' dim3 ') = ' num2str(flatfile(i)) ';']);
    end
end
