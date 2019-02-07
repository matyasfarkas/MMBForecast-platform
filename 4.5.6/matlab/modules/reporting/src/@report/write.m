function o = write(o)
%function o = write(o)
% Write Report object
%
% INPUTS
%   o     [report]  report object
%
% OUTPUTS
%   o     [report]  report object
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

[fid, msg] = fopen(o.fileName, 'w');
if fid == -1
    error(['@report.write: ' msg]);
end

fprintf(fid, '%% Report Object\n');
fprintf(fid, '\\documentclass[11pt]{article}\n');

fprintf(fid, '\\usepackage[%spaper,margin=%f%s', o.paper, o.margin, o.marginUnit);
if strcmpi(o.orientation, 'landscape')
    fprintf(fid, ',landscape');
end
fprintf(fid, ']{geometry}\n');
fprintf(fid, '\\usepackage{pdflscape, booktabs, pgfplots, colortbl, adjustbox, multicol}\n');
fprintf(fid, '\\pgfplotsset{compat=1.5.1}');
fprintf(fid, ['\\makeatletter\n' ...
              '\\def\\blfootnote{\\gdef\\@thefnmark{}\\@footnotetext}\n' ...
              '\\makeatother\n']);

if isoctave && isempty(regexpi(computer, '.*apple.*', 'once'))
    fprintf(fid, '\\usepackage[T1]{fontenc}\n');
    fprintf(fid, '\\usepackage[utf8x]{inputenc}\n');
end
if ispc || ismac
    fprintf(fid, '\\usepgfplotslibrary{fillbetween}\n');
end
fprintf(fid, '\\definecolor{LightCyan}{rgb}{0.88,1,1}\n');
fprintf(fid, '\\definecolor{Gray}{gray}{0.9}\n');
if o.showDate
    fprintf(fid, '\\usepackage{fancyhdr, datetime}\n');
    fprintf(fid, '\\newdateformat{reportdate}{\\THEDAY\\ \\shortmonthname\\ \\THEYEAR}\n');
    fprintf(fid, '\\pagestyle{fancy}\n');
    fprintf(fid, '\\renewcommand{\\headrulewidth}{0pt}\n');
    fprintf(fid, '\\renewcommand{\\footrulewidth}{0.5pt}\n');
    fprintf(fid, '\\rfoot{\\scriptsize\\reportdate\\today\\ -- \\currenttime}\n');
end

% May not need these.....
fprintf(fid, '\\renewcommand{\\textfraction}{0.05}\n');
fprintf(fid, '\\renewcommand{\\topfraction}{0.8}\n');
fprintf(fid, '\\renewcommand{\\bottomfraction}{0.8}\n');
fprintf(fid, '\\setlength{\\parindent}{0in}\n');
fprintf(fid, '\\setlength{\\tabcolsep}{1em}\n');
fprintf(fid, '\\newlength\\sectionheight\n');
if ~isempty(o.header)
    fprintf(fid, '%s\n', o.header);
end
fprintf(fid, '\\begin{document}\n');
if isunix && ~ismac
    fprintf(fid, '\\pgfdeclarelayer{axis background}\n');
    fprintf(fid, '\\pgfdeclarelayer{axis lines}\n');
    fprintf(fid, '\\pgfsetlayers{axis background,axis lines,main}\n');
end
fprintf(fid, '\\pgfplotsset{tick scale binop={\\times},\ntrim axis left}\n');
fprintf(fid, '\\centering\n');

nps = length(o.pages);
for i=1:nps
    if o.showOutput
        fprintf(1, 'Writing Page: %d\n', i);
    end
    o.pages{i}.write(fid, i);
end

fprintf(fid, '\\end{document}\n');
fprintf(fid, '%% End Report Object\n');
status = fclose(fid);
if status == -1
    error('@report.write: closing %s\n', o.fileName);
end
if o.showOutput
    disp('Finished Writing Report!');
end
end
