function [options_, oo_]=ms_compute_mdd(M_, options_, oo_)
% function [options_, oo_]=ms_compute_mdd(M_, options_, oo_)
% Markov-switching SBVAR: Compute Marginal Data Density
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

disp('MS-SBVAR Marginal Data Density');
options_ = set_file_tags(options_);
clean_ms_mdd_files(options_.ms.output_file_tag, options_.ms.proposal_type);
[options_, oo_] = set_ms_estimation_file(options_.ms.file_tag, options_, oo_);
options_ = set_ms_simulation_file(options_);

% setup command line options
opt = ['-mdd -seed ' num2str(options_.DynareRandomStreams.seed)];
opt = [opt ' -ft ' options_.ms.simulation_file_tag];
opt = [opt ' -fto ' options_.ms.output_file_tag];
opt = [opt ' -pf ' options_.ms.mh_file];
opt = [opt ' -d ' num2str(options_.ms.proposal_draws)];
opt = [opt ' -pt ' num2str(options_.ms.proposal_type)];
opt = [opt ' -l '  num2str(options_.ms.proposal_lower_bound)];
opt = [opt ' -h '  num2str(options_.ms.proposal_upper_bound)];
if options_.ms.use_mean_center
    opt = [opt ' -use_mean'];
end

% compute mdd
[err] = ms_sbvar_command_line(opt);
mexErrCheck('ms_compute_mdd',err);

% grab the muller/bridge log mdd from the output file
mull_exp = 'Muller \w+\(\w+\) \= (\d+.\d+e\+\d+)';
bridge_exp = 'Bridge \w+\(\w+\) \= (\d+.\d+e\+\d+)';
bridge_mdd = -1; muller_mdd = -1;
mdd_filename = ['mdd_t' num2str(options_.ms.proposal_type) '_' options_.ms.output_file_tag '.out'];
if exist(mdd_filename,'file')
    mdd_fid = fopen(mdd_filename);
    tline = fgetl(mdd_fid);
    while ischar(tline)
        mull_tok = regexp(tline,mull_exp,'tokens');
        bridge_tok = regexp(tline,bridge_exp,'tokens');
        if (~isempty(mull_tok))
            muller_mdd = str2double(mull_tok{1}{1});
        end
        if (~isempty(bridge_tok))
            bridge_mdd = str2double(bridge_tok{1}{1});
        end
        tline = fgetl(mdd_fid);
    end
    oo_.ms.mueller_log_mdd = muller_mdd;
    oo_.ms.bridged_log_mdd = bridge_mdd;
    fclose(mdd_fid);
end
end
