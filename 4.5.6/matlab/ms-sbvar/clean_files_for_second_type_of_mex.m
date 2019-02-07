function clean_files_for_second_type_of_mex(M_, options_, type)
% function clean_files_for_second_type_of_mex(M_, options_, type)
% clean the files for the appropriate file tag and mex function
%
% INPUTS
%    M_:          (struct)    model structure
%    options_:    (struct)    options
%    type:        (string)    one of irf, forecast or
%    variance_decomposition
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2011-2017 Dynare Team
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

tagtouse = options_.ms.file_tag;
if ~strcmp(tagtouse, options_.ms.output_file_tag)
    tagtouse = options_.ms.output_file_tag;
end

switch type
  case 'irf'
    clean_ms_irf_files(tagtouse);
  case 'forecast'
    clean_ms_forecast_files(tagtouse);
  case 'variance_decomposition'
    clean_ms_variance_decomposition_files(tagtouse);
  otherwise
    error('clean_files_for_second_type_of_mex: should not arrive here');
end
