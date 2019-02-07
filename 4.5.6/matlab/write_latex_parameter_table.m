function write_latex_parameter_table
%function write_latex_parameter_table
% Writes a latex file containing the parameter names, parameter values and
% long names/descriptions
%
% INPUTS
%    none
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

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

global M_

if ~isequal(M_.param_names,M_.param_names_long)
    Long_names_present=1;
else
    Long_names_present=0;
end
fid = fopen([M_.fname '_latex_parameters.tex'], 'w');
fprintf(fid, '\\begin{center}\n');
if Long_names_present==1
    fprintf(fid, '\\begin{longtable}{ccc}\n');
else
    fprintf(fid, '\\begin{longtable}{cc}\n');
end
fprintf(fid, ['\\caption{Parameter Values}\\\\%%\n']);

fprintf(fid, '\\toprule%%\n');
fprintf(fid, '\\multicolumn{1}{c}{\\textbf{Parameter}} &\n');
fprintf(fid, '\\multicolumn{1}{c}{\\textbf{Value}} ');
if Long_names_present==1
    fprintf(fid, '&\n \\multicolumn{1}{c}{\\textbf{Description}}\\\\%%\n');
else
    fprintf(fid, ' \\\\%%\n');
end
fprintf(fid, '\\midrule%%\n');
fprintf(fid, '\\endfirsthead\n');

if Long_names_present==1
    fprintf(fid, '\\multicolumn{3}{c}{{\\tablename} \\thetable{} -- Continued}\\\\%%\n');
else
    fprintf(fid, '\\multicolumn{2}{c}{{\\tablename} \\thetable{} -- Continued}\\\\%%\n');
end
fprintf(fid, '\\midrule%%\n');
fprintf(fid, '\\multicolumn{1}{c}{\\textbf{Parameter}} &\n');
fprintf(fid, '\\multicolumn{1}{c}{\\textbf{Value}} ');
if Long_names_present==1
    fprintf(fid, '&\n  \\multicolumn{1}{c}{\\textbf{Description}}\\\\%%\n');
else
    fprintf(fid, '\\\\%%\n');
end
fprintf(fid, '\\midrule%%\n');
fprintf(fid, '\\endhead\n');

tex = M_.param_names_tex;
long = M_.param_names_long;
for j=1:size(tex,1)
    if Long_names_present==1
        % replace underscores
        long_names_temp=regexprep(strtrim(long(j,:)), '_', '\\_');
        % replace percent
        long_names_temp=regexprep(long_names_temp, '%', '\\%');
        fprintf(fid, '$%s$ \t & \t %4.3f \t & \t %s\\\\\n', ...
                strtrim(tex(j,:)), ...
                M_.params(j,:),...
                long_names_temp);
    else
        fprintf(fid, '$%s$ \t & \t %4.3f \\\\\n', ...
                strtrim(tex(j,:)), ...
                M_.params(j,:));
    end
end
fprintf(fid, '\\bottomrule%%\n');
fprintf(fid, '\\end{longtable}\n');
fprintf(fid, '\\end{center}\n');

fclose(fid);
