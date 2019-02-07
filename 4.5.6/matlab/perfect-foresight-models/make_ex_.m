function oo_=make_ex_(M_,options_,oo_)
% forms oo_.exo_simul and oo_.exo_det_simul
%
% INPUTS
%   M_:           Dynare model structure
%   options_:     Dynare options structure
%   oo_:          Dynare results structure
%
% OUTPUTS
%   oo_:          Dynare results structure
%
% ALGORITHM
%
% SPECIAL REQUIREMENTS
%

% Copyright (C) 1996-2017 Dynare Team
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

global ex0_

if isempty(oo_.exo_steady_state)
    oo_.exo_steady_state = zeros(M_.exo_nbr,1);
end
if M_.exo_det_nbr > 1 && isempty(oo_.exo_det_steady_state)
    oo_.exo_det_steady_state = zeros(M_.exo_det_nbr,1);
end

% Initialize oo_.exo_simul
if isempty(M_.exo_histval)
    if isempty(ex0_)
        oo_.exo_simul = repmat(oo_.exo_steady_state',M_.maximum_lag+options_.periods+M_.maximum_lead,1);
    else
        oo_.exo_simul = [ repmat(ex0_',M_.maximum_lag,1) ; repmat(oo_.exo_steady_state',options_.periods+M_.maximum_lead,1) ];
    end
else
    if isempty(ex0_)
        oo_.exo_simul = [M_.exo_histval'; repmat(oo_.exo_steady_state',options_.periods+M_.maximum_lead,1)];
    else
        error('histval and endval cannot be used simultaneously')
    end
end

% Initialize oo_.exo_det_simul
if M_.exo_det_nbr > 0
    if isempty(M_.exo_det_histval)
        oo_.exo_det_simul = repmat(oo_.exo_det_steady_state',M_.maximum_lag+options_.periods+M_.maximum_lead,1);
    else
        oo_.exo_det_simul = [M_.exo_det_histval'; repmat(oo_.exo_det_steady_state',options_.periods+M_.maximum_lead,1)];
    end
end

% Add temporary shocks
if isfield(M_, 'det_shocks')
    for i = 1:length(M_.det_shocks)
        k = M_.det_shocks(i).periods + M_.maximum_lag;
        ivar = M_.det_shocks(i).exo_id;
        v = M_.det_shocks(i).value;
        if ~M_.det_shocks(i).exo_det
            if ~M_.det_shocks(i).multiplicative
                oo_.exo_simul(k,ivar) = v;
            else
                oo_.exo_simul(k,ivar) = oo_.exo_simul(k,ivar) * v;
            end
        else
            if ~M_.det_shocks(i).multiplicative
                oo_.exo_det_simul(k,ivar) = v;
            else
                oo_.exo_det_simul(k,ivar) = oo_.exo_det_simul(k,ivar) * v;
            end
        end
    end
end
