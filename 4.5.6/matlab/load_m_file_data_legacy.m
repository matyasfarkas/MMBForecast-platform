function o2WysrOISH  = load_m_file_data_legacy(datafile, U7ORsJ0vy3) % --*-- Unitary tests --*--

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

cXDHdrXnqo5KwwVpTRuc6OprAW = datafile(1:end-2);
[pathtocXDHdrXnqo5KwwVpTRuc6OprAW,cXDHdrXnqo5KwwVpTRuc6OprAW,junk] = fileparts(cXDHdrXnqo5KwwVpTRuc6OprAW);

if ~isempty(pathtocXDHdrXnqo5KwwVpTRuc6OprAW)
    % We need to change directory, first we keep the current directory in memory...
    OvMuQsJgjwzYG5Pni0TzU8Acb2YBJva = pwd();
    % Then we move in the directory where the data file is saved.
    cd(pathtocXDHdrXnqo5KwwVpTRuc6OprAW);
end

% We evaluate the matlab script defining the data. All the variables in the
% variables defined in this script are loaded in the current workspace.
eval(cXDHdrXnqo5KwwVpTRuc6OprAW);

if ~isempty(pathtocXDHdrXnqo5KwwVpTRuc6OprAW)
    % If we previously changed directory, we go back to the initial directory.
    cd(OvMuQsJgjwzYG5Pni0TzU8Acb2YBJva);
    clear OvMuQsJgjwzYG5Pni0TzU8Acb2YBJva;
end

% Clear all the variables except the ones defined in the script.
clear('pathtocXDHdrXnqo5KwwVpTRuc6OprAW', 'cXDHdrXnqo5KwwVpTRuc6OprAW', 'junk');

% Get the list of variables in the script.
mj6F4eU1BN = whos();
Z3s1ZFBffw = {mj6F4eU1BN(:).name};

% Check that the variables in varobs are available.
if ~isequal(sort(intersect(Z3s1ZFBffw, U7ORsJ0vy3)), sort(U7ORsJ0vy3))
    qtvUkxKk6b = setdiff(U7ORsJ0vy3, intersect(Z3s1ZFBffw, U7ORsJ0vy3));
    qtvUkxKk6b = sprintf('%s ', qtvUkxKk6b{:});
    qtvUkxKk6b = qtvUkxKk6b(1:end-1);
    error('Some variables are missing (%s)!', qtvUkxKk6b)
end

% Check that the variables are provided as vectors.
N5L9sgRHIu = {};
for uAiwEPcc3Q=1:length(U7ORsJ0vy3)
    if ~isvector(eval(U7ORsJ0vy3{uAiwEPcc3Q}))
        N5L9sgRHIu(end+1) = {U7ORsJ0vy3{uAiwEPcc3Q}};
    end
end
if ~isempty(N5L9sgRHIu)
    N5L9sgRHIu = sprintf('%s ', N5L9sgRHIu{:});
    N5L9sgRHIu = N5L9sgRHIu(1:end-1);
    error('Observed variables should be provided as vectors (%s are not vectors)!')
end

% Check that all the vectors have the same number of elements.
RXZzmKFPFK = numel(eval(U7ORsJ0vy3{1}));
for uAiwEPcc3Q=2:length(U7ORsJ0vy3)
    if ~isequal(numel(eval(U7ORsJ0vy3{1})), RXZzmKFPFK)
        error('All vectors must have the same number of elements (%s has %i elements while %s has %i elements)!', U7ORsJ0vy3{1}, numel(eval(U7ORsJ0vy3{1})), U7ORsJ0vy3{uAiwEPcc3Q}, numel(eval(U7ORsJ0vy3{uAiwEPcc3Q})));
    end
end

% Put the observed variables in data
JSmvfqTSXT = repmat(' vec(%s) ', 1, length(U7ORsJ0vy3));
VbO4y7zOlh = sprintf('[%s]', JSmvfqTSXT);
o2WysrOISH = dseries(eval(sprintf(VbO4y7zOlh, U7ORsJ0vy3{:})), [], U7ORsJ0vy3);

return

%@test:1
% Write a data file
fid = fopen('example1.m','w');
fwriten(fid, 'a = randn(100,1);');
fwriten(fid, 'b = randn(100,1);');
fwriten(fid, 'c = transpose(randn(100,1));');
fwriten(fid, 'd = randn(100,1);');
fwriten(fid, 'e = randn(100,2);');
fwriten(fid, ' ');
fwriten(fid, 'f = NaN(100,1);');
fwriten(fid, 'for i=1:100');
fwriten(fid, '  f(i) = log(rand());')
fwriten(fid, 'end');
fclose(fid);
% Define a set of variables to be loaded.
listofvariablestobeloaded = {'b', 'a'};
% Test if we can load the data.
try
    data = load_m_file_data_legacy('example1.m', listofvariablestobeloaded);
    delete('example1.m');
    t(1) = 1;
catch
    t(1) = 0;
end
T = all(t);
%@eof:1

%@test:2
% Write a data file
fid = fopen('example2.m','w');
fwriten(fid, 'a = randn(100,1);');
fwriten(fid, 'b = randn(100,1);');
fwriten(fid, 'c = transpose(randn(100,1));');
fwriten(fid, 'd = randn(100,1);');
fwriten(fid, 'e = randn(100,2);');
fwriten(fid, ' ');
fwriten(fid, 'f = NaN(100,1);');
fwriten(fid, 'for i=1:100');
fwriten(fid, '  f(i) = log(rand());')
fwriten(fid, 'end');
fclose(fid);
% Define a set of variables to be loaded.
listofvariablestobeloaded = {'e', 'a'};
% Test if we can load the data.
try
    data = load_m_file_data_legacy('example2.m', listofvariablestobeloaded);
    delete('example2.m');
    t(1) = 0;
catch
    t(1) = 1;
end
T = all(t);
%@eof:2

%@test:3
% Write a data file
fid = fopen('example3.m','w');
fwriten(fid, 'a = randn(100,1);');
fwriten(fid, 'b = randn(100,1);');
fwriten(fid, 'c = transpose(randn(100,1));');
fwriten(fid, 'd = randn(100,1);');
fwriten(fid, 'e = randn(100,2);');
fwriten(fid, ' ');
fwriten(fid, 'f = NaN(100,1);');
fwriten(fid, 'for i=1:100');
fwriten(fid, '  f(i) = log(rand());')
fwriten(fid, 'end');
fclose(fid);
% Define a set of variables to be loaded.
listofvariablestobeloaded = {'c', 'a'};
% Test if we can load the data.
try
    data = load_m_file_data_legacy('example3.m', listofvariablestobeloaded);
    delete('example3.m');
    t(1) = 1;
catch
    t(1) = 0;
end
T = all(t);
%@eof:3
