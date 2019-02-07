function o = report_table(varargin)
%function o = report_table(varargin)
% Report_Table Class Constructor
%
% INPUTS
%   0 args => empty report_table
%   1 arg (report_table class) => copy object
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2013-2017 Dynare Team
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

o = struct;

o.tableDirName = 'tmpRepDir';
o.tableName = '';

o.series = {};

o.title = {''};
titleFormatDefalut = {'\large'};
o.titleFormat = titleFormatDefalut;

o.showHlines = false;
o.showVlines = false;
o.vlineAfter = '';
o.vlineAfterEndOfPeriod = false;

o.data = '';
o.seriesToUse = '';
o.range = {};
o.precision = 1;
o.writeCSV = false;

o.highlightRows = {''};

if nargin == 1
    assert(isa(varargin{1}, 'report_table'),['With one arg to Report_Table constructor, ' ...
                        'you must pass a report_table object']);
    o = varargin{1};
    return;
elseif nargin > 1
    if round(nargin/2) ~= nargin/2
        error(['Options to Report_Table constructor must be supplied in name/value ' ...
               'pairs.']);
    end

    optNames = fieldnames(o);

    % overwrite default values
    for pair = reshape(varargin, 2, [])
        ind = find(strcmpi(optNames, pair{1}));
        assert(isempty(ind) || length(ind) == 1);
        if ~isempty(ind)
            o.(optNames{ind}) = pair{2};
        else
            error('%s is not a recognized option to the Report_Table constructor.', pair{1});
        end
    end
end
if ~iscell(o.range)
    o.range = {o.range};
end

if isdates(o.vlineAfter)
    o.vlineAfter = {o.vlineAfter};
end

% Check options provided by user
if ischar(o.title)
    o.title = {o.title};
end
if ischar(o.titleFormat)
    o.titleFormat = {o.titleFormat};
end
if length(o.title) ~= length(o.titleFormat)
    o.titleFormat = repmat(titleFormatDefalut, 1, length(o.title));
end
assert(islogical(o.showHlines), '@report_table.report_table: showHlines must be true or false');
assert(islogical(o.showVlines), '@report_table.report_table: showVlines must be true or false');
assert(isint(o.precision) && o.precision >= 0, '@report_table.report_table: precision must be a non-negative integer');
assert(isempty(o.range) || length(o.range) <=2 && allCellsAreDatesRange(o.range), ...
       ['@report_table.report_table: range is specified as a dates range, e.g. ' ...
        '''dates(''1999q1''):dates(''1999q3'')''.']);
assert(isempty(o.data) || isdseries(o.data), ...
       '@report_table.report_table: data must be a dseries');
assert(isempty(o.seriesToUse) || iscellstr(o.seriesToUse), ...
       '@report_table.report_table: seriesToUse must be a cell array of string(s)');
assert(isempty(o.vlineAfter) || allCellsAreDates(o.vlineAfter), ...
       '@report_table.report_table: vlineAfter must be a dates');
if o.showVlines
    o.vlineAfter = '';
end
assert(islogical(o.vlineAfterEndOfPeriod), ...
       '@report_table.report_table: vlineAfterEndOfPeriod must be true or false');
assert(iscellstr(o.title), ...
       '@report_table.report_table: title must be a cell array of string(s)');
assert(iscellstr(o.titleFormat), ...
       '@report_table.report_table: titleFormat must be a cell array of string(s)');
assert(ischar(o.tableName), '@report_table.report_table: tableName must be a string');
assert(ischar(o.tableDirName), '@report_table.report_table: tableDirName must be a string');
assert(islogical(o.writeCSV), '@report_table.report_table: writeCSV must be either true or false');
assert(iscellstr(o.highlightRows), '@report_table.report_table: highlightRowsmust be a cell string');

% using o.seriesToUse, create series objects and put them in o.series
if ~isempty(o.data)
    if isempty(o.seriesToUse)
        for i=1:o.data.vobs
            o.series{end+1} = report_series('data', o.data{o.data.name{i}});
        end
    else
        for i=1:length(o.seriesToUse)
            o.series{end+1} = report_series('data', o.data{o.seriesToUse{i}});
        end
    end
end
o = rmfield(o, 'seriesToUse');
o = rmfield(o, 'data');

if ~exist(o.tableDirName, 'file')
    mkdir(o.tableDirName);
end

% Create report_table object
o = class(o, 'report_table');
end