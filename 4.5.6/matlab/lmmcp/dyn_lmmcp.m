function [endo_simul,info] = dyn_lmmcp(M,options,oo)

% Copyright (C) 2014 Dynare Team
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

[lb,ub,eq_index] = get_complementarity_conditions(M);

lead_lag_incidence = M.lead_lag_incidence;

ny = M.endo_nbr;

max_lag = M.maximum_endo_lag;

nyp = nnz(lead_lag_incidence(1,:)) ;
iyp = find(lead_lag_incidence(1,:)>0) ;
ny0 = nnz(lead_lag_incidence(2,:)) ;
iy0 = find(lead_lag_incidence(2,:)>0) ;
nyf = nnz(lead_lag_incidence(3,:)) ;
iyf = find(lead_lag_incidence(3,:)>0) ;

nd = nyp+ny0+nyf;
nrc = nyf+1 ;
isp = [1:nyp] ;
is = [nyp+1:ny+nyp] ;
isf = iyf+nyp ;
isf1 = [nyp+ny+1:nyf+nyp+ny+1] ;
stop = 0 ;
iz = [1:ny+nyp+nyf];

periods = options.periods;
steady_state = oo.steady_state;
params = M.params;
endo_simul = oo.endo_simul;
exo_simul = oo.exo_simul;
i_cols_1 = nonzeros(lead_lag_incidence(2:3,:)');
i_cols_A1 = find(lead_lag_incidence(2:3,:)');
i_cols_T = nonzeros(lead_lag_incidence(1:2,:)');
i_cols_j = 1:nd;
i_upd = ny+(1:periods*ny);

x = endo_simul(:);

model_dynamic = str2func([M.fname,'_dynamic']);
z = x(find(lead_lag_incidence'));
[res,A] = model_dynamic(z, exo_simul, params, steady_state,2);
nnzA = nnz(A);

LB = repmat(lb,periods,1);
UB = repmat(ub,periods,1);

Y0 = endo_simul(:,1);
YT = endo_simul(:,end);
x = endo_simul(:,2:end-1);
x = x(:);

func_handle = @(x) dyn_lmmcp_func(x,model_dynamic, Y0, YT, exo_simul, ...
                                  params, steady_state, periods, ny, ...
                                  lead_lag_incidence, i_cols_A1, i_cols_1, ...
                                  i_cols_T, i_cols_j,nnzA,eq_index);

[x, info] = lmmcp(func_handle,x,LB,UB,options.lmmcp);

endo_simul = [Y0 reshape(x,ny,periods) YT];