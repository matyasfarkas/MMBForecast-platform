% Script to set the necessary paths for the Reporting toolbox

% Copyright (C) 2015-2017 Dynare Team
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

% Find reporting source directory
reporting_src_root = [strrep(which('initialize_reporting_toolbox'),'initialize_reporting_toolbox.m','') 'src'];

% Add path to reporting source
addpath(reporting_src_root);
addpath([reporting_src_root filesep '..' filesep 'macros']);

% Reminder to add and initialize dates & dseries toolboxes
if ~exist('emptydatesobject', 'var')
    disp('Remember to add the paths to the dates toolbox before working with the reporting toolbox');
end

if ~exist('emptydseriesobject', 'var')
    disp('Remember to add the paths to the dseries toolbox before working with the reporting toolbox');
end
