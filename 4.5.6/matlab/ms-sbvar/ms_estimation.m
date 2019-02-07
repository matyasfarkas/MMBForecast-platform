function [options_, oo_]=ms_estimation(M_, options_, oo_)
% function [options_, oo_]=ms_estimation(M_, options_, oo_)
% Markov-switching SBVAR: Estimation
%
% INPUTS
%    M_:          (struct)    model structure
%    options_:    (struct)    options
%    oo_:         (struct)    results
%
% OUTPUTS
%    options_:    (struct)    options
%    oo_:         (struct)    results
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

disp('MS-SBVAR Estimation');
options_ = set_file_tags(options_);

% general setup
if options_.ms.create_init
    clean_ms_init_files(options_.ms.file_tag);
    ms_sbvar_setup(options_);
    clean_ms_estimation_files(options_.ms.file_tag);
    clean_ms_estimation_files(options_.ms.output_file_tag);
else
    if exist(['init_' options_.ms.file_tag '.dat'],'file') ~= 2 || ...
            exist(['matlab_' options_.ms.file_tag '.prn'],'file') ~= 2 || ...
            exist([options_.ms.file_tag '_markov_file.dat'],'file') ~= 2
        error(['You must run ms_estimation without no_create_init as '...
               'you are missing an initialization file it needs.']);
    end
    if ~strcmp(options_.ms.file_tag, options_.ms.output_file_tag)
        clean_ms_estimation_files(options_.ms.output_file_tag);
    end
end

% setup command line options
opt = ['-estimate -seed ' num2str(options_.DynareRandomStreams.seed)];
opt = [opt ' -ft ' options_.ms.file_tag];
opt = [opt ' -fto ' options_.ms.output_file_tag];
opt = [opt ' -cb ' num2str(options_.ms.convergence_starting_value)];
opt = [opt ' -ce ' num2str(options_.ms.convergence_ending_value)];
opt = [opt ' -ci ' num2str(options_.ms.convergence_increment_value)];
opt = [opt ' -ib ' num2str(options_.ms.max_iterations_starting_value)];
opt = [opt ' -ii ' num2str(options_.ms.max_iterations_increment_value)];
opt = [opt ' -mb ' num2str(options_.ms.max_block_iterations)];
opt = [opt ' -repeat_max ' num2str(options_.ms.max_repeated_optimization_runs)];
opt = [opt ' -repeat_tol_obj ' num2str(options_.ms.function_convergence_criterion)];
opt = [opt ' -repeat_tol_parms ' num2str(options_.ms.parameter_convergence_criterion)];
opt = [opt ' -random ' num2str(options_.ms.number_of_large_perturbations)];
opt = [opt ' -random_small ' num2str(options_.ms.number_of_small_perturbations)];
opt = [opt ' -random_small_ndraws ' num2str(options_.ms.number_of_posterior_draws_after_perturbation)];
opt = [opt ' -random_max ' num2str(options_.ms.max_number_of_stages)];
opt = [opt ' -random_tol_obj ' num2str(options_.ms.random_function_convergence_criterion)];
opt = [opt ' -random_tol_parms ' num2str(options_.ms.random_parameter_convergence_criterion)];

% estimation
[err] = ms_sbvar_command_line(opt);
mexErrCheck('ms_estimation', err);

[options_, oo_] = set_ms_estimation_file(options_.ms.output_file_tag, options_, oo_);
end
