function plan = basic_plan(plan, exogenous, expectation_type, date, value)
% Adds a simple shock to the forecast scenario plan
%
% INPUTS
%  o plan                 [structure]       A structure describing the different shocks and the implied variables, the date of the shocks and the path of the shock (forecast scenario).
%                                           The plan structure is created by the functions init_plan, basic_plan and flip_plan
%  o exogenous            [string]          A string containing the name of the exognous shock.
%  o expectation_type     [string]          A string indicating the type of expectation: 'surprise' for an unexpected shock, and 'perfect_foresight' for a perfectly anticpated shock.
%  o date                 [dates]           The period of the shock
%  o value                [array of double] A vector of double containing the values of the exogenous variable
%
%
% OUTPUTS
%  plan                   [structure]        Returns a structure containing the updated forecast scenario.
%
%
% Copyright (C) 2013-2017 Dynare Team
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

if ~ischar(expectation_type) || size(expectation_type,1) ~= 1
    error(['in basic_plan the third argument should be a string containing the simulation type (''perfect_foresight'' or ''surprise'')']);
end
exogenous = strtrim(exogenous);
ix = find(strcmp(exogenous, plan.exo_names));
if  isempty(ix)
    error(['in basic_plan the second argument ' exogenous ' is not an exogenous variable']);
end
sdate = length(date);
if sdate > 1
    if date(1) < plan.date(1) || date(end) > plan.date(end)
        error(['in basic_plan the fourth argument (date='  date ') must lay inside the plan.date ' plan.date]);
    end
else
    if date < plan.date(1) || date > plan.date(end)
        error(['in basic_plan the fourth argument (date='  date ') must lay iside the plan.date ' plan.date]);
    end
end
if length(date) ~= length(value)
    error(['in basic_plan the number of dates (' int2str(length(date)) ') is not equal to the numbers of shock (' int2str(length(value)) ') for exogenous variable ' exogenous]);
end
if ~isempty(plan.options_cond_fcst_.controlled_varexo)
    common_var = find(ix == plan.options_cond_fcst_.controlled_varexo);
    if ~isempty(common_var)
        common_date = intersect(date, plan.constrained_date_{common_var});
        if ~isempty(common_date)
            [date_, i_date] = setdiff(date, common_date);
            value = value(i_date);
            if common_date.length > 1
                the_dates = [cell2mat(strings(common_date(1))) ':' cell2mat(strings(common_date(end)))];
            else
                the_dates = cell2mat(strings(common_date));
            end
            warning(['Impossible case: ' plan.exo_names{plan.options_cond_fcst_.controlled_varexo(common_var)} ' is used both as a shock and as an endogenous variable to control the path of ' plan.endo_names{plan.constrained_vars_(common_var)} ' at the dates ' the_dates]);
            warning('This shock will not be considered');
        end
    end
end
if isempty(plan.shock_vars_)
    plan.shock_vars_ = ix;
    if strcmp(expectation_type, 'perfect_foresight')
        plan.shock_perfect_foresight_ = 1;
    else
        plan.shock_perfect_foresight_ = 0;
    end
else
    plan.shock_vars_ = [plan.shock_vars_ ; ix];
    if strcmp(expectation_type, 'perfect_foresight')
        plan.shock_perfect_foresight_ = [plan.shock_perfect_foresight_ ; 1];
    else
        plan.shock_perfect_foresight_ = [plan.shock_perfect_foresight_ ; 0];
    end
end
plan.shock_date_{length(plan.shock_date_) + 1} = date;
plan.shock_str_date_{length(plan.shock_str_date_) + 1} = strings(date);
plan.shock_int_date_{length(plan.shock_int_date_) + 1} = date - plan.date(1) + 1;
plan.shock_paths_{length(plan.shock_paths_) + 1} = value;
