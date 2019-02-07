function dr=set_state_space(dr,DynareModel,DynareOptions)
% Write the state space representation of the reduced form solution.

%@info:
%! @deftypefn {Function File} {[@var{dr} =} set_state_space (@var{dr},@var{DynareModel},@var{DynareOptions})
%! @anchor{set_state_space}
%! @sp 1
%! Write the state space representation of the reduced form solution.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item dr
%! Matlab's structure describing decision and transition rules.
%! @item DynareModel
%! Matlab's structure describing the model (initialized by dynare, see @ref{M_})
%! @item DynareOptions
%! Matlab's structure describing the current options (initialized by dynare, see @ref{options_}).
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dr
%! Matlab's structure describing decision and transition rules.
%! @end table
%! @sp 2
%! @strong{This function is called by:}
%! @sp 1
%! @ref{check}, @ref{discretionary_policy_1}, @ref{dynare_estimation_init}, @ref{dyn_risky_steady_state_solver}, @ref{osr1}, @ref{partial_information/dr1_PI}, @ref{pea/pea_initialization}, @ref{stochastic_solvers}, @ref{stoch_simul}
%! @sp 2
%! @strong{This function calls:}
%! @sp 2
%! @end deftypefn
%@eod:

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

max_lead = DynareModel.maximum_endo_lead;
max_lag = DynareModel.maximum_endo_lag;
endo_nbr = DynareModel.endo_nbr;
lead_lag_incidence = DynareModel.lead_lag_incidence;
klen = max_lag + max_lead + 1;

fwrd_var = find(lead_lag_incidence(max_lag+2:end,:))';
if max_lag > 0
    pred_var = find(lead_lag_incidence(1,:))';
    both_var = intersect(pred_var,fwrd_var);
    pred_var = setdiff(pred_var,both_var);
    fwrd_var = setdiff(fwrd_var,both_var);
    stat_var = setdiff([1:endo_nbr]',union(union(pred_var,both_var),fwrd_var));  % static variables
else
    pred_var = [];
    both_var = [];
    stat_var = setdiff([1:endo_nbr]',fwrd_var);
end
if DynareOptions.block == 1
    order_var = DynareModel.block_structure.variable_reordered;
else
    order_var = [ stat_var(:); pred_var(:); both_var(:); fwrd_var(:)];
end
inv_order_var(order_var) = (1:endo_nbr);

% building kmask for z state vector in t+1
if max_lag > 0
    kmask = [];
    if max_lead > 0
        kmask = lead_lag_incidence(max_lag+2,order_var) ;
    end
    kmask = [kmask; lead_lag_incidence(1,order_var)] ;
else
    if max_lead==0 %%in this case lead_lag_incidence has no entry max_lag+2
        error('Dynare currently does not allow to solve purely static models in a stochastic context.')
    end
    kmask = lead_lag_incidence(max_lag+2,order_var) ;
end

kmask = kmask';
kmask = kmask(:);
i_kmask = find(kmask);
nd = nnz(kmask);           % size of the state vector
kmask(i_kmask) = (1:nd);
% auxiliary equations

% composition of state vector
% col 1: variable;           col 2: lead/lag in z(t+1);
% col 3: A cols for t+1 (D); col 4: A cols for t (E)
kstate = [ repmat([1:endo_nbr]',klen-1,1) kron([klen:-1:2]',ones(endo_nbr,1)) ...
           zeros((klen-1)*endo_nbr,2)];
kiy = flipud(lead_lag_incidence(:,order_var))';
kiy = kiy(:);
if max_lead > 0
    kstate(1:endo_nbr,3) = kiy(1:endo_nbr)-nnz(lead_lag_incidence(max_lag+1,:));
    kstate(kstate(:,3) < 0,3) = 0;
    kstate(endo_nbr+1:end,4) = kiy(2*endo_nbr+1:end);
else
    kstate(:,4) = kiy(endo_nbr+1:end);
end
kstate = kstate(i_kmask,:);

dr.order_var = order_var;
dr.inv_order_var = inv_order_var';
dr.kstate = kstate;

dr.transition_auxiliary_variables = [];
