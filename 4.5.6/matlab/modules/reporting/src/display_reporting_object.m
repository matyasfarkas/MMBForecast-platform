function display_reporting_object(o)
%function display_reporting_object(o)
% Display a Reporting Object
%   i.e., one of: graph
%                 page
%                 report
%                 report_series
%                 report_table
%                 section
%                 vspace
%
% INPUTS
%   o   [object] reporting object
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

fprintf('\n%s Object = \n\n', upper(class(o)));
fields = fieldnames(o);
for i=1:length(fields)
    fprintf('    %s: ', fields{i});
    val = o.(fields{i});
    if iscell(val)
        fprintf('{');
        for j=1:length(val)
            if ischar(val{j});
                fprintf('''%s''', val{j});
                if j~=length(val)
                    fprintf(', ');
                end
            else
                if strcmp(fields{i}, 'pages') || strcmp(fields{i}, 'sections')
                    fprintf('%d', length(val));
                    break;
                else
                    fprintf('fix this');
                end
            end
        end
        fprintf('}');
    elseif ischar(val)
        fprintf('''%s''', val);
    elseif isnumeric(val)
        fprintf('%s', num2str(val));
    elseif islogical(val)
        if val
            fprintf('true');
        else
            fprintf('false');
        end
    elseif isobject(val)
        if isdates(val)
            if isempty(val)
                fprintf('<dates: empty>');
            else
                fprintf('<dates: %s, ..., %s>', ...
                        date2string(val(1)), date2string(val(end)));
            end
        elseif isdseries(val)
            if numel(val) == 1
                fprintf('<dseries: %s>', val.name{1});
            else
                fprintf('%s', class(val));
            end
        else
            cl = class(val);
            fprintf('%d', val.(['num' upper(cl(1)) cl(2:end)]));
        end
    else
        fprintf('fix this');
    end
    fprintf('\n');
end
fprintf('\n');
end