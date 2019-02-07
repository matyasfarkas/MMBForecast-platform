function [options_, oo_]=ms_simulation(M_, options_, oo_)
% function [options_, oo_]=ms_simulation(M_, options_, oo_)
% Markov-switching SBVAR: Simulation
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

disp('MS-SBVAR Simulation');
options_ = set_file_tags(options_);
clean_ms_simulation_files(options_.ms.output_file_tag);
[options_, oo_] = set_ms_estimation_file(options_.ms.file_tag, options_, oo_);

% setup command line options
opt = ['-simulate -seed ' num2str(options_.DynareRandomStreams.seed)];
opt = [opt ' -ft ' options_.ms.file_tag];
opt = [opt ' -fto ' options_.ms.output_file_tag];
opt = [opt ' -ndraws ' num2str(options_.ms.mh_replic)];
opt = [opt ' -burnin ' num2str(options_.ms.drop)];
opt = [opt ' -thin ' num2str(options_.ms.thinning_factor)];
opt = [opt ' -mh ' num2str(options_.ms.adaptive_mh_draws)];

if options_.ms.save_draws
    opt = [opt ' -flat '];
end

% simulation
[err] = ms_sbvar_command_line(opt);
mexErrCheck('ms_simulation',err);
end
