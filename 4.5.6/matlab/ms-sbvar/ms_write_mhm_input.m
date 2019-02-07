function ms_write_mhm_input(fname, options_ms)
% function ms_write_mhm_input(fname, options_ms)
%
% INPUTS
%    fname:         (string)    filename
%    options_ms:    (struct)    options
%
% OUTPUTS
%    none
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

fh = fopen(fname,'w');


fprintf(fh,'/**********************************************************\n');
fprintf(fh,' *** This input file is read by swzmsbvar_mhm_1 and swzmsbvar_mhm_1.exe only, NOT by swzmsbvar_printdraws.exe.\n');
fprintf(fh,' ***\n');
fprintf(fh,' **********************************************************/\n');

fprintf(fh,'\n\n//------------- 1st set of posterior draws to find optimal scales for Metropolis (30000). ---------------\n');
fprintf(fh,'//== number draws for first burn-in ==//  //For determining the Metropolis scales only.\n');
fprintf(fh,'%d\n\n',options_ms.draws_nbr_burn_in_1);

fprintf(fh,'//------------- MCMC burn-in draws once the Metropolis scales (previous stage) are fixed. --------------\n');
fprintf(fh,'//------------- 2nd set of standard burn-in posterior draws to throw away the initial draws (10000). ---------------\n');
fprintf(fh,'//== number draws for second burn-in ==//\n');
fprintf(fh,'%d\n\n',options_ms.draws_nbr_burn_in_2);

fprintf(fh,'//--------------- 1st set of posterior draws to compute the mean and variance for the weighting function in the MHM (200000) ----------------\n');
fprintf(fh,'//== number draws to estimate mean and variance ==//\n');
fprintf(fh,'%d\n\n',options_ms.draws_nbr_mean_var_estimate);

fprintf(fh,'//--------------- Only applied to mhm_2 process: total number of MCMC draws = thinning factor * 2nd set of saved posterior draws ----------------\n');
fprintf(fh,'//== thinning factor for modified harmonic mean process ==//\n');
fprintf(fh,'%d\n\n',options_ms.thinning_factor);

fprintf(fh,'//--------------- 2nd set of saved posterior draws from MHM_2 (second stage): saved draws AFTER thinning (1000000) ----------------\n');
fprintf(fh,'//== number draws for modified harmonic mean process ==//\n');
fprintf(fh,'%d\n\n',options_ms.draws_nbr_modified_harmonic_mean);

fprintf(fh,'//------- 1st stage: computing all three tightness factors for Dirichlet.  ---------\n');
fprintf(fh,'//------- 2nd stage: hard-code the second scale factor (in principle, we can do all three). ---------\n');
fprintf(fh,'//------- It seems that Dan''s code only use the first element of the following scales.  The scale applies to the Dirichlet''s hyperparameter alpha for the diagonal of the transition matrix in the weighting function.  Note that the weighting function for the transition matrix parameters is Dirichlet. ---------\n');

fprintf(fh,'//== scale values for Dirichlet distribution ==//\n');
fprintf(fh,'3\n\n');
fprintf(fh,'%f ',options_ms.dirichlet_scale);
fprintf(fh,'\n');
fclose(fh);