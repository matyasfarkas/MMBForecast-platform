function histvalf(fname)

% Copyright (C) 2014-2017 Dynare Team
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

global M_ oo_ ex0_

if ~exist(fname)
    error(['Can''t find datafile: ' fname ]);
end

M_.endo_histval = repmat(oo_.steady_state, 1, M_.maximum_endo_lag);

% Also fill in oo_.exo_simul: necessary if we are in deterministic context,
% since aux vars for lagged exo are not created in this case
if isempty(oo_.exo_simul)
    if isempty(ex0_)
        oo_.exo_simul = repmat(oo_.exo_steady_state',M_.maximum_lag,1);
    else
        oo_.exo_simul = repmat(ex0_',M_.maximum_lag,1);
    end
end

S = load(fname);

outvars = fieldnames(S);

for i = 1:length(outvars)
    ov_ = outvars{i};
    if ov_(end) == '_'
        ov = ov_(1:end-1);
        j = strmatch(ov, M_.endo_names, 'exact');
        if isempty(j)
            warning(['smoother2histval: output variable ' ov ' does not exist.'])
        end
    else
        % Lagged endogenous or exogenous, search through aux vars
        undidx = find(ov_ == '_', 1, 'last'); % Index of last underscore in name
        ov = ov_(1:(undidx-1));
        lead_lag = ov_((undidx+1):end);
        lead_lag = regexprep(lead_lag,'l','-');
        lead_lag = str2num(lead_lag);
        j = [];
        for i = 1:length(M_.aux_vars)
            if M_.aux_vars(i).type ~= 1 && M_.aux_vars(i).type ~= 3
                continue
            end
            if M_.aux_vars(i).type == 1
                % Endogenous
                orig_var = deblank(M_.endo_names(M_.aux_vars(i).orig_index, :));
            else
                % Exogenous
                orig_var = deblank(M_.exo_names(M_.aux_vars(i).orig_index, :));
            end
            if strcmp(orig_var, ov) && M_.aux_vars(i).orig_lead_lag == lead_lag
                j = M_.aux_vars(i).endo_index;
            end
        end
        if isempty(j)
            % There is no aux var corresponding to (orig_var, lead_lag).
            % If this is an exogenous variable, then it means we should put
            % the value in oo_.exo_simul (we are probably in deterministic
            % context).
            k = strmatch(ov, M_.exo_names);
            if isempty(k)
                warning(['smoother2histval: output variable ' ov '(' lead_lag ') does not exist.'])
            else
                oo_.exo_simul((M_.maximum_lag-M_.maximum_endo_lag+1):M_.maximum_lag, k) = getfield(S, ov_);
            end
            continue
        end
    end
    M_.endo_histval(j, :) = getfield(S, ov_);
end
