function [endogenousvariables, exogenousvariables] = model_inversion(constraints, ...
                                                  exogenousvariables, ...
                                                  initialconditions, DynareModel, DynareOptions, DynareOutput)

% INPUTS
% - constraints         [dseries]        with N constrained endogenous variables from t1 to t2.
% - exogenousvariables  [dseries]        with Q exogenous variables.
% - initialconditions   [dseries]        with M endogenous variables starting before t1 (M initialcond must contain at least the state variables).
% - DynareModel         [struct]         M_, Dynare global structure containing informations related to the model.
% - DynareOptions       [struct]         options_, Dynare global structure containing all the options.
%
% OUTPUTS
% - endogenous          [dseries]
% - exogenous           [dseries]
%
% REMARKS

% Copyright (C) 2017 Dynare Team
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

if ~isequal(nargin, 6)
    error('model_inversion: This routine require six input arguments!')
end

if ~isdseries(constraints)
    error('model_inversion: First input argument must be a dseries object!')
end

if ~isdseries(exogenousvariables)
    error('model_inversion: Second input argument must be a dseries object!')
end

if ~isdseries(initialconditions)
    error('model_inversion: Third input argument must be a dseries object!')
end

if ~isstruct(DynareModel)
    error('model_inversion: Last input argument must be structures (M_)!')
end

% Set range where the endogenous variables are constrained.
crange = constraints.dates;

% Check that the number of instruments match the number of constrained endogenous variables.
instruments = exogenousvariables(crange);
freeinnovations = instruments.name(find(all(isnan(instruments))));
if ~isequal(length(freeinnovations), constraints.vobs)
    error('model_inversion: The number of instruments must be equal to the number of constrained variables!')
end

% Check if some of the exogenous variables are given.
observed_exogenous_variables_flag = false;
if exogenousvariables.vobs>constraints.vobs
    observed_exogenous_variables_flag = true;
end

% Get the list of endogenous and exogenous variables.
endo_names = cellstr(DynareModel.endo_names);
exo_names = cellstr(DynareModel.exo_names);

% Use specidalized routine if the model is backward looking.
if ~DynareModel.maximum_lead
    [endogenousvariables, exogenousvariables] = ...
        backward_model_inversion(constraints, exogenousvariables, initialconditions, ...
                                 endo_names, exo_names, freeinnovations, ...
                                 DynareModel, DynareOptions, DynareOutput);
    return
end

% Initialize fplan
fplan = init_plan(crange);

% Set the exogenous observed variables.
if observed_exogenous_variables_flag
    list_of_observed_exogenous_variables = setdiff(exo_names, freeinnovations);
    observed_exogenous_variables = exogenousvariables{list_of_observed_exogenous_variables{:}};
    for i=1:length(list_of_observed_exogenous_variables)
        fplan = basic_plan(fplan, list_of_observed_exogenous_variables{i}, ...
                           'surprise', crange, observed_exogenous_variables{list_of_observed_exogenous_variables{i}}.data(2:length(crange)+1));
    end
end

% Set constrained path for the endogenous variables.
for i = 1:constraints.vobs
    fplan = flip_plan(fplan, constraints.name{i}, freeinnovations{i}, 'surprise', crange, transpose(constraints.data(:,i)));
end

% Identify the innovations (model inversion)
f = det_cond_forecast(fplan, initialconditions, crange);

endogenousvariables = f{endo_names{:}};
exogenousvariables = f{exo_names{:}};