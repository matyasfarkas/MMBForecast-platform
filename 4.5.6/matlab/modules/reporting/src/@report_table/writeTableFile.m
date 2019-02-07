function o = writeTableFile(o, pg, sec, row, col)
%function o = writeTableFile(o, pg, sec, row, col)
% Write a Report_Table object
%
% INPUTS
%   o   [report_table]    report_table object
%   pg  [integer] this page number
%   sec [integer] this section number
%   row [integer] this row number
%   col [integer] this col number
%
% OUTPUTS
%   o   [report_table]    report_table object
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

ne = length(o.series);
if ne == 0
    warning('@report_table.write: no series to plot, returning');
    return;
end

if isempty(o.tableName)
    o.tableName = sprintf('%s/table_pg%d_sec%d_row%d_col%d.tex', o.tableDirName, pg, sec, row, col);
else
    o.tableName = [o.tableDirName '/' o.tableName];
end

[fid, msg] = fopen(o.tableName, 'w');
if fid == -1
    error(['@report_table.writeTableFile: ' msg]);
end

%number of left-hand columns, 1 until we allow the user to group data,
% e.g.: GDP Europe
%         GDP France
%         GDP Germany
% this example would be two lh columns, with GDP Europe spanning both
nlhc = 1;

if isempty(o.range)
    dates = getMaxRange(o.series);
    o.range = {dates};
else
    dates = o.range{1};
end
ndates = dates.ndat;

fprintf(fid, '%% Report_Table Object\n');
fprintf(fid, '\\setlength{\\parindent}{6pt}\n');
fprintf(fid, '\\setlength{\\tabcolsep}{4pt}\n');
fprintf(fid, '\\begin{tabular}{@{}l');

for i=1:ndates
    fprintf(fid, 'r');
    if o.showVlines
        fprintf(fid, '|');
    elseif o.vlineAfterEndOfPeriod && dates(i).time(2) == dates(i).freq
        fprintf(fid, '|');
    elseif ~isempty(o.vlineAfter)
        for j=1:length(o.vlineAfter)
            if dates(i) == o.vlineAfter{j}
                fprintf(fid, '|');
            end
        end
    end
end
datedata = dates.time;
years = unique(datedata(:, 1));
if length(o.range) > 1
    rhscols = strings(o.range{2});
    if o.range{2}.freq == 1
        rhscols = strrep(rhscols, 'Y', '');
    end
else
    rhscols = {};
end
for i=1:length(rhscols)
    fprintf(fid, 'r');
    if o.showVlines
        fprintf(fid, '|');
    end
end
nrhc = length(rhscols);
ncols = ndates+nlhc+nrhc;
fprintf(fid, '@{}}%%\n');
for i=1:length(o.title)
    if ~isempty(o.title{i})
        fprintf(fid, '\\multicolumn{%d}{c}{%s %s}\\\\\n', ...
                ncols, o.titleFormat{i}, o.title{i});
    end
end
fprintf(fid, '\\toprule%%\n');

% Column Headers
thdr = num2cell(years, size(years, 1));
if dates.freq == 1
    for i=1:size(thdr, 1)
        fprintf(fid, ' & %d', thdr{i, 1});
    end
    for i=1:length(rhscols)
        fprintf(fid, ' & %s', rhscols{i});
    end
else
    thdr{1, 2} = datedata(:, 2)';
    if size(thdr, 1) > 1
        for i=2:size(thdr, 1)
            split = find(thdr{i-1, 2} == dates.freq, 1, 'first');
            assert(~isempty(split), '@report_table.writeTableFile: Shouldn''t arrive here');
            thdr{i, 2} = thdr{i-1, 2}(split+1:end);
            thdr{i-1, 2} = thdr{i-1, 2}(1:split);
        end
    end
    for i=1:size(thdr, 1)
        fprintf(fid, ' & \\multicolumn{%d}{c}{%d}', size(thdr{i,2}, 2), thdr{i,1});
    end
    for i=1:length(rhscols)
        fprintf(fid, ' & %s', rhscols{i});
    end
    fprintf(fid, '\\\\\n');
    switch dates.freq
      case 4
        sep = 'Q';
      case 12
        sep = 'M';
      case 52
        sep = 'W';
      otherwise
        error('@report_table.writeTableFile: Invalid frequency.');
    end
    for i=1:size(thdr, 1)
        period = thdr{i, 2};
        for j=1:size(period, 2)
            fprintf(fid, ' & \\multicolumn{1}{c');
            if o.showVlines
                fprintf(fid, '|');
            elseif o.vlineAfterEndOfPeriod && j == size(period, 2)
                fprintf(fid, '|');
            elseif ~isempty(o.vlineAfter)
                for k=1:length(o.vlineAfter)
                    if o.vlineAfter{k}.time(1) == thdr{i} && ...
                            o.vlineAfter{k}.time(2) == period(j)
                        fprintf(fid, '|');
                    end
                end
            end
            fprintf(fid, '}{%s%d}', sep, period(j));
        end
    end
end
fprintf(fid, '\\\\[-2pt]%%\n');
fprintf(fid, '\\hline%%\n');
fprintf(fid, '%%\n');

% Write Report_Table Data
if o.writeCSV
    csvseries = dseries();
end
for i=1:ne
    o.series{i}.writeSeriesForTable(fid, o.range, o.precision, ncols, o.highlightRows{mod(i,length(o.highlightRows))+1});
    if o.writeCSV
        if isempty(o.series{i}.tableSubSectionHeader)
            csvseries = [csvseries ...
                         o.series{i}.data(dates).set_names([...
                             num2str(i) '_' ...
                             o.series{i}.data.name{:}])];
        end
    end
    if o.showHlines
        fprintf(fid, '\\hline\n');
    end
end
if o.writeCSV
    csvseries.save(strrep(o.tableName, '.tex', ''), 'csv');
end
fprintf(fid, '\\bottomrule\n');
fprintf(fid, '\\end{tabular}\\setlength{\\parindent}{0pt}\n \\par \\medskip\n\n');
fprintf(fid, '%% End Report_Table Object\n');
if fclose(fid) == -1
    error('@report_table.writeTableFile: closing %s\n', o.filename);
end
end
