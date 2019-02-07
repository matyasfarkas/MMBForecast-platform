function ms_write_markov_file(fname, options)
% function ms_write_markov_file(fname, options)
%
% INPUTS
%    fname:       (string)    name of markov file
%    options:     (struct)    options
%
% OUTPUTS
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

n_chains = length(options.ms.ms_chain);
nvars = length(options.varobs);

fh = fopen(fname,'w');
%/******************************************************************************/
%/********************* Markov State Variable Information **********************/
%/******************************************************************************/

fprintf(fh,'//== Markov State Variables with Simple Restrictions ==//\n\n');


%//This number is NOT used but read in.
fprintf(fh,'//== Number Observations ==//\n');
fprintf(fh,'0\n\n');

fprintf(fh,'//== Number independent state variables ==//\n');
fprintf(fh,'%d\n\n',n_chains);

for i_chain = 1:n_chains

    %//=====================================================//
    %//== state_variable[i] (1 <= i <= n_state_variables) ==//
    %//=====================================================//
    fprintf(fh,'//== Number of states for state_variable[%d] ==//\n', ...
            i_chain);
    n_states = length(options.ms.ms_chain(i_chain).regime);
    fprintf(fh,'%d\n\n',n_states);

    %//== 03/15/06: DW TVBVAR code reads the data below and overwrite the prior data read somewhere else if any.
    %//== Each column contains the parameters for a Dirichlet prior on the corresponding
    %//== column of the transition matrix.  Each element must be positive.  For each column,
    %//== the relative size of the prior elements determine the relative size of the elements
    %//== of the transition matrix and overall larger sizes implies a tighter prior.
    fprintf(fh,['//== Transition matrix prior for state_variable[%d] ==//\n'], ...
            i_chain);
    Alpha = ones(n_states,n_states);
    for i_state = 1:n_states
        p = 1-1/options.ms.ms_chain(i_chain).regime(i_state).duration;
        Alpha(i_state,i_state) = p*(n_states-1)/(1-p);
        fprintf(fh,'%22.16f',Alpha(i_state,:));
        fprintf(fh,'\n');
    end

    fprintf(fh,['\n//== Dirichlet dimensions for state_variable[%d] ' ...
                '==//\n'],i_chain);
    %        fprintf(fh,'%d ',repmat(n_states,1,n_states));
    fprintf(fh,'%d ',repmat(2,1,n_states));
    fprintf(fh,'\n\n');

    %//== The jth restriction matrix is n_states-by-free[j].  Each row of the restriction
    %//== matrix has exactly one non-zero entry and the sum of each column must be one.
    fprintf(fh,['//== Column restrictions for state_variable[%d] ' ...
                '==//\n'],i_chain);
    for i_state = 1:n_states
        if i_state == 1
            M = eye(n_states,2);
        elseif i_state == n_states
            M = [zeros(n_states-2,2); eye(2)];
        else
            M = zeros(n_states,2);
            M(i_state+[-1 1],1) = ones(2,1)/2;
            M(i_state,2) = 1;
            disp(M)
        end
        for j_state = 1:n_states
            fprintf(fh,'%f ',M(j_state,:));
            fprintf(fh,'\n');
        end
        fprintf(fh,'\n');
    end
end

%/******************************************************************************/
%/******************************* VAR Parameters *******************************/
%/******************************************************************************/
%//NOT read
fprintf(fh,'//== Number Variables ==//\n');
fprintf(fh,'%d\n\n',nvars);

%//NOT read
fprintf(fh,'//== Number Lags ==//\n');
fprintf(fh,'%d\n\n',options.ms.nlags);

%//NOT read
fprintf(fh,'//== Exogenous Variables ==//\n');
fprintf(fh,'1\n\n');


%//== nvar x n_state_variables matrix.  In the jth row, a non-zero value implies that
%this state variable controls the jth column of A0 and Aplus
fprintf(fh,['//== Controlling states variables for coefficients ==//\' ...
            'n']);

for i_var = 1:nvars
    for i_chain = 1:n_chains
        if ~isfield(options.ms.ms_chain(i_chain),'svar_coefficients') ...
                || isempty(options.ms.ms_chain(i_chain).svar_coefficients)
            i_equations = 0;
        else
            i_equations = ...
                options.ms.ms_chain(i_chain).svar_coefficients.equations;
        end
        if strcmp(i_equations,'ALL') || any(i_equations == i_var)
            fprintf(fh,'%d ',1);
        else
            fprintf(fh,'%d ',0);
        end
    end
    fprintf(fh,'\n');
end

%//== nvar x n_state_variables matrix.  In the jth row, a non-zero value implies that
%this state variable controls the jth diagonal element of Xi
fprintf(fh,'\n//== Controlling states variables for variance ==//\n');
for i_var = 1:nvars
    for i_chain = 1:n_chains
        if ~isfield(options.ms.ms_chain(i_chain),'svar_variances') ...
                || isempty(options.ms.ms_chain(i_chain).svar_variances)
            i_equations = 0;
        else
            i_equations = ...
                options.ms.ms_chain(i_chain).svar_variances.equations;
        end
        if strcmp(i_equations,'ALL') || any(i_equations == i_var)
            fprintf(fh,'%d ',1);
        else
            fprintf(fh,'%d ',0);
        end
    end
    fprintf(fh,'\n');
end

fclose(fh);
