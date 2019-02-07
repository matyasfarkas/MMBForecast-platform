function dyn_saveas(h, fname, nodisplay, graph_format)
%function dyn_saveas(h, fname, nodisplay, graph_format)
% save figures for DYNARE
%
% INPUTS
%    h     : figure handle
%    fname : name of the saved figure
%    nodisplay: the value of the command-specific nodisplay argument or options_.nodisplay
%    graph_format: the value of the command-specific graph_format argument or options_.graph_format
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2012-2017 Dynare Team
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

if any(strcmp('eps',cellstr(graph_format)))
    if isoctave
        fname = strrep(fname,'/',filesep);
        fname = strrep(fname,'\',filesep);
        if nodisplay && ispc
            set(h, 'Visible','on');
        end
    end
    print(h,'-depsc2',[fname,'.eps'])
end
if any(strcmp('pdf',cellstr(graph_format)))
    if isoctave
        error('Octave cannot create pdf files!')
    else
        print(h,'-dpdf',[fname,'.pdf'])
    end
end
if any(strcmp('fig',cellstr(graph_format)))
    if isoctave
        error('Octave cannot create fig files!')
    else
        if nodisplay
            %  THE FOLLOWING LINES COULD BE USED IF BUGS/PROBLEMS ARE REPORTED USING LINE 60
            %             set(h,'Units','Normalized')
            %             mypos=get(h,'Position');
            %             set(h,'Position',[-1 -1 mypos(3:4)])
            %             set(h, 'Visible','on');
            set(h,'CreateFcn','set(gcf, ''Visible'',''on'')') ;
        end
        saveas(h,[fname '.fig']);
    end
end
if any(strcmp('none',cellstr(graph_format)))
    % don't save
    % check here as a reminder that none is an option to graph_format
end
if nodisplay
    close(h);
end
