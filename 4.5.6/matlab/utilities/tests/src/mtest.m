function [check, info] = mtest(fname, fpath)

% Extracts unit test sections from matlab's routine, executes the tests and reports results.
%
% INPUTS
%  - fname [string], name of the Matlab routine where unit tests have to be run.
%  * fpath [string], path to the routine
%
% OUTPUTS
%  - check [integer], scalar equal to 0 if the test fails and 0 otherwise
%  - info  [cell], a cell describing the test results. Cell info has nn rows and
%          five columns. Each row correponds to a unitary test in fname, and the
%          columns report the following informations:
%
%            Column 1 Name of the tested routine.
%            Column 2 Number of the unitary test.
%            Column 3 Status of the unitary test (0 if the unitary test fails, 1 otherwise).
%            Column 4 Details about the failure (vector of 0 and 1).
%            Column 5 Elapsed time in seconds.
%
% REMARKS
%  - If only one input argument is provided, fname must be a string containing the
%    full path to the targeted matlab routine.
%  - If two input arguments are provided, fname is the base name of the targeted
%    matlab routine and fpath is the path to this routine.

% Copyright (C) 2013-2017 Dynare Team
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

% Default answer (no problem).
check = 1;

% Open the matlab file.
if nargin<2 || isempty(fpath)
    fname = append_extension_if_needed(fname);
    fid = fopen(fname,'r');
    [junk, FNAME, vessel] = fileparts(fname);
else
    fname = append_extension_if_needed(fname);
    fid = fopen([fpath filesep fname],'r');
    FNAME = fname(1:end-2);
end

% Read the matlab file.
file = textscan(fid,'%s','delimiter','\n');
file = file{1};

% Close the matlab file.
fclose(fid);

% Locate the test blocks.
b1 = find(strncmp(file,'%@test:',7))+1;
b2 = find(strncmp(file,'%@eof:',6))-1;
nn = length(b2);

if length(b1)-length(b2)
    error('test:: There is a problem with the test blocks definition!')
end

% Initialize the second output if necessary.
if nargout>1
    info = cell(nn,5);
end

% Perform the tests.
for i=1:nn
    if nargout>1
        info(i,1) = {fname};
        info(i,2) = {i};
    end
    % Write the temporary test routine.
    tid = fopen([FNAME '_test_' int2str(i) '.m'],'w');
    fprintf(tid,['function [T,t,LOG] = ' FNAME '_test_' int2str(i) '()\n']);
    fprintf(tid,['try\n']);
    if (length(file{b1(i)+1})>2 && isequal(file{b1(i)+1}(1:3), '%$ ')) || (length(file{b1(i)+1})>1 && isequal(file{b1(i)+1}(1:2), '%$'))
        remove_first_columns = true;
    else
        remove_first_columns = false;
    end
    for j=b1(i):b2(i)
        if remove_first_columns
            str = sprintf('%s \n',file{j}(4:end));
        else
            str = sprintf('%s \n',file{j}(1:end));
        end
        str = regexprep(str, '%', '%%');
        fprintf(tid,str);
    end
    fprintf(tid,['LOG = NaN;\n']);
    if isoctave
        fprintf(tid,'catch\n');
        fprintf(tid,'exception = lasterror;\n');
        fprintf(tid, 'LOG = ''%s'';\n','The Log output is not available with Octave!');
    else
        fprintf(tid,'catch exception\n');
        fprintf(tid,['LOG = getReport(exception,''extended'');\n']);
    end
    fprintf(tid,['T = NaN;\n']);
    fprintf(tid,['t = NaN;\n']);
    fprintf(tid,['end\n']);
    fclose(tid);
    % Call the temporary test routine.
    tic;
    [TestFlag,TestDetails,LOG] = feval([FNAME '_test_' int2str(i)]);
    time = toc;
    if isnan(TestFlag)
        fprintf(['\n'])
        fprintf(['Call to ' FNAME ' test routine n°' int2str(i) ' failed (' datestr(now) ')!\n'])
        fprintf(['\n'])
        if ~isoctave
            disp(LOG)
        end
        check = 0;
        if nargout>1
            info(i,3) = {0};
        end
        continue
    end
    if ~TestFlag
        if nargout>1
            info(i,3) = {0};
            tmp = ones(length(TestDetails),1);
        end
        fprintf(['Test n°' int2str(i) ' for routine ' FNAME ' failed (' datestr(now) ')!\n']);
        for j=1:length(TestDetails)
            if ~TestDetails(j)
                if nargout>1
                    tmp(j) = 0;
                end
                fprintf(['Output argument n°' int2str(j) ' didn''t give the expected result.\n']);
            end
        end
        if nargout>1
            info(i,4) = {tmp};
            info(i,5) = {NaN};
        end
        check = 0;
    else
        if nargout>1
            info(i,3) = {1};
            info(i,4) = {ones(length(TestDetails),1)};
            info(i,5) = {time};
        end
        delete([FNAME '_test_' int2str(i) '.m'])
    end
end


function file = append_extension_if_needed(file)
if ~isequal(file(end-1:end),'.m')
    file = [file '.m'];
end