function clean_ms_mdd_files(file_tag, pt)
% function clean_ms_mdd_files(file_tag, pt)
% removes MS mdd files
%
% INPUTS
%    file_tag: string indicating tag to use when deleting files
%    pt:       numeric proposal type, used to create filename
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2011 Dynare Team
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

delete_if_exists(['mdd_t' num2str(pt) '_' file_tag '.out']);
delete_if_exists(['proposal_t' num2str(pt) '_' file_tag '.out']);
end
