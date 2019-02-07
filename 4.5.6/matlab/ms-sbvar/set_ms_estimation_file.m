function [options_, oo_]=set_ms_estimation_file(file_tag, options_, oo_)
% function [options_, oo_]=set_ms_estimation_file(file_tag, options_, oo_)
% Set oo_.ms.maxparams, oo_.ms.A0, oo_.ms.Aplus, oo_.ms.Zeta, oo_.ms.Q
% based on estimation output files
%
% INPUTS
%    file_tag:    (string)    necessary because of different meanings of
%                             file_tag between ms_estimation and other ms_*
%                             routines
%    options_:    (struct)    options
%    oo_:         (struct)    results
%
% OUTPUTS
%    options_:    (struct)    options
%    oo_:         (struct)    results
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2011-2012 Dynare Team
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

options_.ms.free_param_file = ['est_free_' file_tag '.out'];
if ~exist(options_.ms.free_param_file,'file')
    error(['ERROR: Could not find free parameter file: ' options_.ms.free_param_file]);
else
    oo_.ms.maxparams = load(options_.ms.free_param_file);
    oo_.ms.maxparams = oo_.ms.maxparams(3:end)';
end

options_.ms.VAR_parameters_file = [file_tag '.mat'];
if ~exist(options_.ms.VAR_parameters_file,'file')
    error(['ERROR: Could not find VAR parameters file: ' options_.ms.VAR_parameters_file]);
else
    oo_.ms = load(options_.ms.VAR_parameters_file);
end
end
