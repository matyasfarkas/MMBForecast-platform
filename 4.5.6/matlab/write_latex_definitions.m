function write_latex_definitions
%function write_latex_definitions
% Writes a latex file containing the variable names, latex names, and
% tags/comments
%
% INPUTS
%    none
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

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

global M_;

if M_.exo_det_nbr == 0
    tables = {'Endogenous', 'Exogenous', 'Parameters'};
    M_var_root = {'M_.endo', 'M_.exo', 'M_.param'};
else
    tables = {'Endogenous', 'Exogenous', 'Exogenous Deterministic', 'Parameters'};
    M_var_root = {'M_.endo', 'M_.exo', 'M_.exo_det', 'M_.param'};
end

fid = fopen([M_.fname '_latex_definitions.tex'], 'w');
for i=1:length(tables)
    fprintf(fid, '\\begin{center}\n');
    fprintf(fid, '\\begin{longtable}{ccc}\n');
    fprintf(fid, ['\\caption{' tables{i} '}\\\\%%\n']);

    fprintf(fid, '\\hline%%\n');
    fprintf(fid, '\\multicolumn{1}{c}{\\textbf{Variable}} &\n');
    fprintf(fid, '\\multicolumn{1}{c}{\\textbf{\\LaTeX}} &\n');
    fprintf(fid, '\\multicolumn{1}{c}{\\textbf{Description}}\\\\%%\n');
    fprintf(fid, '\\hline\\hline%%\n');
    fprintf(fid, '\\endfirsthead\n');

    fprintf(fid, '\\multicolumn{3}{c}{{\\tablename} \\thetable{} -- Continued}\\\\%%\n');
    fprintf(fid, '\\hline%%\n');
    fprintf(fid, '\\multicolumn{1}{c}{\\textbf{Variable}} &\n');
    fprintf(fid, '\\multicolumn{1}{c}{\\textbf{\\LaTeX}} &\n');
    fprintf(fid, '\\multicolumn{1}{c}{\\textbf{Description}}\\\\%%\n');
    fprintf(fid, '\\hline\\hline%%\n');
    fprintf(fid, '\\endhead\n');

    names = eval([M_var_root{i} '_names']);
    tex = eval([M_var_root{i} '_names_tex']);
    long = eval([M_var_root{i} '_names_long']);
    for j=1:size(names,1)
        fprintf(fid, '\\texttt{%s} & $%s$ & %s\\\\\n', ...
                regexprep(strtrim(names(j,:)), '_', '\\_'), ...
                strtrim(tex(j,:)), ...
                regexprep(strtrim(long(j,:)), '_', '\\_'));
    end
    fprintf(fid, '\\hline%%\n');
    fprintf(fid, '\\end{longtable}\n');
    fprintf(fid, '\\end{center}\n');
end
fclose(fid);
end
