function o = write(o, fid, pg)
%function o = write(o, fid, pg)
% Write a Page object
%
% INPUTS
%   o              [page]     page object
%   fid            [integer]  file id
%   pg             [integer]  this page number
%
% OUTPUTS
%   o              [page]     page object
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

assert(fid ~= -1);

fprintf(fid, '\n%% Page Object\n');
if strcmpi(o.orientation, 'landscape')
    fprintf(fid, '\\begin{landscape}\n');
end

for i=1:length(o.footnote)
    fprintf(fid, '\\blfootnote{\\tiny %d. %s}', i, o.footnote{i});
end
fprintf(fid,'\n');

if ~isempty(o.latex)
    if ~exist(o.pageDirName, 'dir')
        mkdir(o.pageDirName)
    end
    pagename = [o.pageDirName '/page_' num2str(pg) '.tex'];
    [fidp, msg] = fopen(pagename, 'w');
    if fidp == -1
        error(['@page.write: ' msg]);
    end
    fprintf(fidp, '%s', o.latex);
    if fclose(fidp) == -1
        error('@page.write: closing %s\n', pagename);
    end
    fprintf(fid, '\\input{%s}', pagename);
else
    fprintf(fid, '\\begin{tabular}[t]{c}\n');
    for i=1:length(o.title)
        if isint(o.titleTruncate)
            if length(o.title{i}) > o.titleTruncate
                o.title{i} = o.title{i}(1:o.titleTruncate);
            end
        end
        fprintf(fid,'\\multicolumn{1}{c}{%s %s}\\\\\n', o.titleFormat{i}, o.title{i});
    end

    nps = length(o.sections);
    for i=1:nps
        o.sections{i}.write(fid, pg, i);
    end
    fprintf(fid, '\\end{tabular}\n');
end

if strcmpi(o.orientation, 'landscape')
    fprintf(fid, '\\end{landscape}\n');
end
fprintf(fid, '\\clearpage\n');
fprintf(fid, '%% End Page Object\n\n');
end
