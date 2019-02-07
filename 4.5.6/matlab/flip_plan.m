function plan = flip_plan(plan, exogenous, endogenous, expectation_type, date, value)
% Adds to the forecast scenario a conditional forecast shock (the path of an endogenous variable is constrained and the values compatible values of the related exogenous variable will be compued)
%
% INPUTS
%  o plan                 [structure]       A structure describing the different shocks and the implied variables, the date of the shocks and the path of the shock (forecast scenario).
%                                           The plan structure is created by the functions init_plan, basic_plan and flip_plan
%  o exogenous            [string]          A string containing the name of the endogenous variable with a constrained path.
%  o endogenous           [string]          A string containing the name of the exogenous. This exogenous variable will be en endogenous variable when the conditional forecast will be perform.
%  o expectation_type     [string]          A string indicating the type of expectation: 'surprise' for an unexpected shock, and 'perfect_foresight' for a perfectly anticpated shock.
%  o date                 [dates]           The period of the shock
%  o value                [array of double] A vector of double containing the constrined path on the endogenous variable
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
    error(['in flip_plan the fourth argument should be a string containing the simulation type (''perfect_foresight'' or ''surprise'')']);
end
exogenous = strtrim(exogenous);
ix = find(strcmp(exogenous, plan.endo_names));
if  isempty(ix)
    error(['in flip_plan the second argument ' exogenous ' is not an endogenous variable']);
end
endogenous = strtrim(endogenous);
iy = find(strcmp(endogenous, plan.exo_names));
if  isempty(iy)
    error(['in flip_plan the third argument ' endogenous ' is not an exogenous variable']);
end
sdate = length(date);
if sdate > 1
    if date(1) < plan.date(1) || date(end) > plan.date(end)
        error(['in flip_plan the fifth argument (date='  date ') must lay inside the plan.date ' plan.date]);
    end
else
    if date < plan.date(1) || date > plan.date(end)
        error(['in flip_plan the fifth argument (date='  date ') must lay iside the plan.date ' plan.date]);
    end
end
if ~isempty(plan.shock_vars_)
    common_var = find(iy == plan.shock_vars_);
    if ~isempty(common_var)
        common_date = intersect(date, plan.shock_date_{common_var});
        if ~isempty(common_date)
            if common_date.length > 1
                the_dates = [cell2mat(strings(common_date(1))) ':' cell2mat(strings(common_date(end)))];
            else
                the_dates = cell2mat(strings(common_date));
            end
            error(['Impossible case: ' plan.exo_names{plan.shock_vars_(common_var)} ' is used both as a shock and as an endogenous variable to control the path of ' plan.endo_names{ix} ' at the dates ' the_dates]);
        end
    end
end
i_ix = find(ix == plan.constrained_vars_);
if isempty(i_ix)
    if isempty(plan.constrained_vars_)
        plan.constrained_vars_ = ix;
        plan.options_cond_fcst_.controlled_varexo  = iy;
        if strcmp(expectation_type, 'perfect_foresight')
            plan.constrained_perfect_foresight_ = 1;
        else
            plan.constrained_perfect_foresight_ = 0;
        end
    else
        plan.constrained_vars_ = [plan.constrained_vars_ ; ix];
        plan.options_cond_fcst_.controlled_varexo  = [plan.options_cond_fcst_.controlled_varexo ; iy];
        if strcmp(expectation_type, 'perfect_foresight')
            plan.constrained_perfect_foresight_ = [plan.constrained_perfect_foresight_ ; 1];
        else
            plan.constrained_perfect_foresight_ = [plan.constrained_perfect_foresight_ ; 0];
        end
    end
    plan.constrained_date_{length(plan.constrained_date_) + 1} = date;
    plan.constrained_str_date_{length(plan.constrained_str_date_) + 1} = strings(date);
    plan.constrained_int_date_{length(plan.constrained_int_date_) + 1} = date - plan.date(1) + 1;
    plan.constrained_paths_{length(plan.constrained_paths_) + 1} = value;
elseif plan.options_cond_fcst_.controlled_varexo(i_ix) == iy % same exogenous and endogenous hard tune
[plan.constrained_str_date_{i_ix}, i1, i2] = union(strings(date), plan.constrained_str_date_{i_ix});
plan.constrained_date_{i_ix} = [date(i1) plan.constrained_date_{i_ix}(i2)];
plan.constrained_int_date_{i_ix} = [date(i1) - plan.date(1) + 1; plan.constrained_int_date_{i_ix}(i2)];
plan.constrained_paths_{i_ix} = [value(i1)'; plan.constrained_paths_{i_ix}(i2)];
else
    error(['impossible case you have two conditional forecasts:\n - one involving ' plan.endo_names{plan.options_cond_fcst_.controlled_varexo(i_ix),:} ' as control and ' plan_exo_names{plan.constrained_vars_(ix_)} ' as constrined endogenous\n - the other involving  ' plan.endo_names{plan.options_cond_fcst_.controlled_varexo(iy),:} ' as control and ' plan_exo_names{plan.constrained_vars_(ix)} ' as constrined endogenous\n']);
end
