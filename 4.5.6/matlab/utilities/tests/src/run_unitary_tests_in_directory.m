function report = run_unitary_tests_in_directory(dirname, savereport, printreport, sendreport)

% Runs all the unitary tests defined in a directory (and subfolders).
%
% INPUTS
%  - dirname     [string], name of the directory where
%  - savereport  [integer], scalar equal to 0 or 1. If equal to 1 the report is saved in a mat file.
%  - printreport [integer], scalar equal to 0 or 1. If equal to 1 the report is printed on screen.
%  - sendreport  [integer], scalar equal to 0 or 1. If equal to 1 the report is sent by email.
%
% OUTPUTS
%  - report      [cell], first output argument of run_unitary_test routine.
%
% REMARKS
%  1. Git needs to be available on the system and it is assumed that the content of dirname is versionned with Git.
%
% See also get_directory_description, run_unitary_tests, build_report_summary

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

INIT_PATH = pwd();

cd(dirname);

system('git show --pretty=format:"Last commit %H by %an, %ar %n-> %s" HEAD > git.info');
system('git rev-parse HEAD > git.last-commit-hash');

fid = fopen('git.info');
gitinfo = fgetl(fid);
gitinfo = char(gitinfo,fgetl(fid));
fclose(fid);

fid = fopen('git.last-commit-hash');
gitlastcommithash = fgetl(fid);
fclose(fid);

cd(INIT_PATH);

matlabverion = version;
platform = computer;

listoffiles = get_directory_description(dirname);

diary(['report-' gitlastcommithash '.log'] )

str = sprintf('Unitary tests in %s', dirname);
lstr = length(str);
sstr = repmat('*', 1, lstr);
skipline()
disp(sstr)
disp(str)
disp(sstr)

[report, time] = run_unitary_tests(listoffiles);
diary off

if nargin>1 && savereport>0
    save(['report-' gitlastcommithash '.mat'],'report','time','gitinfo','gitlastcommithash','matlabverion','platform');
end

if nargin>2
    build_report_summary(['report-' gitlastcommithash '.mat'], printreport, sendreport);
end