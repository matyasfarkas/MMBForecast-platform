function results = prior_sampler(drsave,M_,bayestopt_,options_,oo_,estim_params_)
% This function builds a (big) prior sample.
%
% INPUTS
%   drsave      [integer]    Scalar. If equal to 1, then dr structure is saved with each prior draw.
%   M_          [structure]  Model description.
%   bayestopt_  [structure]  Prior distribution description.
%   options_    [structure]  Global options of Dynare.
%
% OUTPUTS:
%   results     [structure]  Various statistics.
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2009-2016 Dynare Team
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

% Initialization.
prior_draw(bayestopt_, options_.prior_trunc);
PriorDirectoryName = CheckPath('prior/draws',M_.dname);
work = ~drsave;
iteration = 0;
loop_indx = 0;
file_indx = [];
count_bk_indeterminacy = 0;
count_bk_unstability = 0;
count_bk_singularity = 0;
count_static_var_def = 0;
count_no_steadystate = 0;
count_steadystate_file_exit = 0;
count_dll_problem = 0;
count_complex_jacobian = 0;
count_complex_steadystate = 0;
count_nan_steadystate = 0;
count_nan_params = 0;
count_complex_params = 0;
count_nonpositive_steadystate = 0;
count_unknown_problem = 0;
NumberOfSimulations = options_.prior_mc;
NumberOfParameters = length(bayestopt_.p1);
NumberOfEndogenousVariables = size(M_.endo_names,1);
NumberOfElementsPerFile = ceil(options_.MaxNumberOfBytes/NumberOfParameters/NumberOfEndogenousVariables/8) ;

if NumberOfSimulations <= NumberOfElementsPerFile
    TableOfInformations = [ 1 ,  NumberOfSimulations , 1] ;
    NumberOfFiles = 1;
else
    NumberOfFiles = ceil(NumberOfSimulations/NumberOfElementsPerFile) ;
    NumberOfElementsInTheLastFile = NumberOfSimulations - NumberOfElementsPerFile*(NumberOfFiles-1) ;
    TableOfInformations = NaN(NumberOfFiles,3) ;
    TableOfInformations(:,1) = transpose(1:NumberOfFiles) ;
    TableOfInformations(1:NumberOfFiles-1,2) = NumberOfElementsPerFile*ones(NumberOfFiles-1,1) ;
    TableOfInformations(NumberOfFiles,2) = NumberOfElementsInTheLastFile ;
    TableOfInformations(1,3) = 1;
    TableOfInformations(2:end,3) = cumsum(TableOfInformations(1:end-1,2))+1;
end

if drsave
    pdraws = cell(TableOfInformations(1,2),drsave+2) ;
else
    pdraws = NaN(TableOfInformations(1,2),NumberOfParameters+1);
end

sampled_prior_expectation = zeros(NumberOfParameters,1);
sampled_prior_covariance  = zeros(NumberOfParameters,NumberOfParameters);

file_line_number = 0;
file_indx_number = 0;

oo_.dr=set_state_space(oo_.dr,M_,options_);

hh = dyn_waitbar(0,'Please wait. Prior sampler...');
set(hh,'Name','Prior sampler.');

% Simulations.
while iteration < NumberOfSimulations
    if ~mod(iteration,10)
        dyn_waitbar(iteration/NumberOfSimulations,hh,'Please wait. Prior sampler...');
    end
    loop_indx = loop_indx+1;
    params = prior_draw();
    M_ = set_all_parameters(params,estim_params_,M_);
    [dr,INFO,M_,options_,oo_] = resol(work,M_,options_,oo_);
    file_line_number = file_line_number + 1;
    iteration = iteration + 1;
    if drsave
        pdraws(file_line_number,1) = {params};
        pdraws(file_line_number,2) = {INFO(1)};
    else
        pdraws(file_line_number,1:NumberOfParameters) = params;
        pdraws(file_line_number,NumberOfParameters+1) = INFO(1);
    end
    switch INFO(1)
      case 0
        if drsave
            pdraws(file_line_number,3) = {dr};
        end
        [sampled_prior_expectation,sampled_prior_covariance] = ...
            recursive_prior_moments(sampled_prior_expectation,sampled_prior_covariance,params,iteration);
      case 1
        count_static_undefined = count_static_undefined + 1;
      case 2
        count_dll_problem = count_dll_problem + 1;
      case 3
        count_bk_unstability = count_bk_unstability + 1 ;
      case 4
        count_bk_indeterminacy = count_bk_indeterminacy + 1 ;
      case 5
        count_bk_singularity = count_bk_singularity + 1 ;
      case 20
        count_no_steadystate = count_no_steadystate + 1 ;
      case 19
        count_steadystate_file_exit = count_steadystate_file_exit + 1 ;
      case 6
        count_complex_jacobian = count_complex_jacobian + 1 ;
      case 21
        count_complex_steadystate = count_complex_steadystate + 1 ;
      case 22
        count_nan_steadystate = count_nan_steadystate + 1 ;
      case 23
        count_complex_params = count_complex_params + 1 ;
      case 24
        count_nan_params = count_nan_params + 1 ;
      case 26
        count_nonpositive_steadystate = count_nonpositive_steadystate + 1;
      otherwise
        count_unknown_problem = count_unknown_problem + 1 ;
    end
    if ( file_line_number==TableOfInformations(file_indx_number+1,2) )
        file_indx_number = file_indx_number + 1;
        save([ PriorDirectoryName '/prior_draws' int2str(file_indx_number) '.mat' ],'pdraws');
        if file_indx_number<NumberOfFiles
            if drsave
                pdraws = cell(TableOfInformations(file_indx_number+1,2),drsave+2);
            else
                pdraws = NaN(TableOfInformations(file_indx_number+1,2),NumberOfParameters+1);
            end
        end
        file_line_number = 0;
    end
end

dyn_waitbar_close(hh);

% Get informations about BK conditions and other things...
results.bk.indeterminacy_share = count_bk_indeterminacy/loop_indx;
results.bk.unstability_share = count_bk_unstability/loop_indx;
results.bk.singularity_share = count_bk_singularity/loop_indx;
results.dll.problem_share = count_dll_problem/loop_indx;
results.ss.problem_share = count_no_steadystate/loop_indx;
results.ss.complex_share = count_complex_steadystate/loop_indx;
results.ass.problem_share = count_steadystate_file_exit/loop_indx;
results.ss.nonpositive_share = count_nonpositive_steadystate/loop_indx;
results.jacobian.problem_share = count_complex_jacobian/loop_indx;
results.garbage_share = ...
    results.bk.indeterminacy_share + ...
    results.bk.unstability_share + ...
    results.bk.singularity_share + ...
    results.dll.problem_share + ...
    results.ss.problem_share + ...
    results.ass.problem_share + ...
    results.ss.complex_share + ...
    results.ss.nonpositive_share + ...
    results.jacobian.problem_share + ...
    count_unknown_problem/loop_indx ;
results.prior.mean = sampled_prior_expectation;
results.prior.variance = sampled_prior_covariance;
results.prior.mass = 1-results.garbage_share;

function [mu,sigma] = recursive_prior_moments(m0,s0,newobs,iter)
%  Recursive estimation of order one and two moments (expectation and
%  covariance matrix). newobs should be a row vector. I do not use the
%  function recursive_moments here, because this function is to be used when
%  newobs is a 2D array.
m1 = m0 + (newobs'-m0)/iter;
qq = m1*m1';
s1 = s0 + ( (newobs'*newobs-qq-s0) + (iter-1)*(m0*m0'-qq') )/iter;
mu = m1;
sigma = s1;