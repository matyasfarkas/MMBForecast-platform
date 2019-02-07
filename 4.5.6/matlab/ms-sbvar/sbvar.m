function sbvar(M, options)
% function sbvar(M, options)
%
% INPUTS
%    M_:          (struct)    model structure
%    options_:    (struct)    options
%
% OUTPUTS
%   none.
%
% ALGORITHM
%   ...
%
% SPECIAL REQUIREMENTS
%   none.
%

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

clean_sbvar_files();
options.data = read_variables(options.datafile,options.varobs,[],options.xls_sheet,options.xls_range);
[options.ms.final_year,options.ms.final_subperiod] = check_datafile_years_assigned(options);

if options.forecast == 0
    options.forecast = 4;
end

if options.ms.upper_cholesky
    if options.ms.lower_cholesky
        error(['Upper Cholesky and lower Cholesky decomposition can''t be ' ...
               'requested at the same time!'])
    else
        options.ms.restriction_fname = 'upper_cholesky';
    end
elseif options.ms.lower_cholesky
    options.ms.restriction_fname = 'lower_cholesky';
elseif ~isempty(options.ms.Qi) && ~isempty(options.ms.Ri)
    options.ms.restriction_fname = 'exclusions';
end

ms_mardd(options);
end
