function dynInfo(fun)

%@info:
%! @deftypefn {Function File} dynInfo (@var{fun})
%! @anchor{dynInfo}
%! @sp 1
%! Displays internal documentation of matlab/octave routine @var{fun}.m.
%! @sp 2
%!
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item fun
%! string, name of the matlab/octave routine for which internal documentation is needed.
%! @end table
%! @sp 2
%!
%! @strong{Outputs}
%! @sp 1
%! None.
%! @sp 2
%!
%! @strong{This function is called by:}
%! @sp 1
%! @ref{internals}, @ref{build_internal_documentation}
%!
%! @strong{This function calls:}
%! @sp 1
%! @ref{get_internal_doc_block}.
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2011-2013 Dynare Team
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

% Original author: stephane DOT adjemian AT univ DASH lemans DOT fr

if isempty(strfind(fun,'@')) & (~isempty(strfind(fun,'/')) || ~isempty(strfind(fun,'\')) )
    [pathstr1, name, ext] = fileparts(fun);
    addpath(pathstr1);
    rm_path = 1;
else
    rm_path = 0;
end

[pathstr2, name, ext] = fileparts(which(fun));

if strcmp(ext(2:end),'m')
    block = get_internal_doc_block(name,pathstr2);
    if ~isempty(block)
        fid = fopen([fun '.texi'],'wt');
        for i=1:size(block,1)
            fprintf(fid,'%s\n',deblank(block(i,:)));
        end
        fclose(fid);
        skipline(2)
        system(['makeinfo --plaintext --no-split --no-validate ' fun '.texi']);
        delete([fun '.texi']);
    else
        disp('No documentation for this routine!')
    end
else
    disp('Not a known matlab/octave routine!')
end

if rm_path
    rmpath(pathstr1)
end