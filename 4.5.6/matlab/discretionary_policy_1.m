function [dr,ys,info]=discretionary_policy_1(oo_,Instruments)

% Copyright (C) 2007-2017 Dynare Team
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

global M_ options_
persistent Hold

dr = [];
ys = [];
info = 0;

if isempty(options_.qz_criterium)
    options_.qz_criterium = 1+1e-6;
end

% safeguard against issues like running ramsey policy first and then running discretion
if isfield(M_,'orig_model')
    orig_model = M_.orig_model;
    M_.endo_nbr = orig_model.endo_nbr;
    M_.endo_names = orig_model.endo_names;
    M_.lead_lag_incidence = orig_model.lead_lag_incidence;
    M_.maximum_lead = orig_model.maximum_lead;
    M_.maximum_endo_lead = orig_model.maximum_endo_lead;
    M_.maximum_lag = orig_model.maximum_lag;
    M_.maximum_endo_lag = orig_model.maximum_endo_lag;
else
    M_.orig_model = M_;
end

beta = get_optimal_policy_discount_factor(M_.params,M_.param_names);

exo_nbr = M_.exo_nbr;
if isfield(M_,'orig_model')
    orig_model = M_.orig_model;
    endo_nbr = orig_model.endo_nbr;
    endo_names = orig_model.endo_names;
    lead_lag_incidence = orig_model.lead_lag_incidence;
    MaxLead = orig_model.maximum_lead;
    MaxLag = orig_model.maximum_lag;
else
    endo_names = M_.endo_names;
    endo_nbr = M_.endo_nbr;
    MaxLag=M_.maximum_lag;
    MaxLead=M_.maximum_lead;
    lead_lag_incidence = M_.lead_lag_incidence;
end

%call steady_state_file if present to update parameters
if options_.steadystate_flag
    % explicit steady state file
    [junk,M_.params,info] = evaluate_steady_state_file(oo_.steady_state,[oo_.exo_steady_state; oo_.exo_det_steady_state],M_, ...
                                                      options_,0);
end
[U,Uy,W] = feval([M_.fname,'_objective_static'],zeros(endo_nbr,1),[], M_.params);
if any(any(Uy~=0))
    error(['discretionary_policy: the objective function must have zero ' ...
           'first order derivatives'])
end

W=reshape(W,endo_nbr,endo_nbr);

klen = MaxLag + MaxLead + 1;
iyv=lead_lag_incidence';
% Find the jacobian
z = repmat(zeros(endo_nbr,1),1,klen);
z = z(nonzeros(iyv)) ;
it_ = MaxLag + 1 ;

if exo_nbr == 0
    oo_.exo_steady_state = [] ;
end

[junk,jacobia_] = feval([M_.fname '_dynamic'],z, [zeros(size(oo_.exo_simul)) ...
                    oo_.exo_det_simul], M_.params, zeros(endo_nbr,1), it_);
if any(junk~=0)
    error(['discretionary_policy: the model must be written in deviation ' ...
           'form and not have constant terms'])
end

eq_nbr= size(jacobia_,1);
instr_nbr=endo_nbr-eq_nbr;

if instr_nbr==0
    error('discretionary_policy:: There are no available instruments, because the model has as many equations as variables.')
end
if size(Instruments,1)< instr_nbr
    error('discretionary_policy:: There are fewer declared instruments than omitted equations.')
elseif size(Instruments,1)> instr_nbr
    error('discretionary_policy:: There are more declared instruments than omitted equations.')    
end

instr_id=nan(instr_nbr,1);
for j=1:instr_nbr
    vj=deblank(Instruments(j,:));
    vj_id=strmatch(vj,endo_names,'exact');
    if ~isempty(vj_id)
        instr_id(j)=vj_id;
    else
        error([mfilename,':: instrument ',vj,' not found'])
    end
end

Indices={'lag','0','lead'};
iter=1;
for j=1:numel(Indices)
    eval(['A',Indices{j},'=zeros(eq_nbr,endo_nbr);'])
    if strcmp(Indices{j},'0')||(strcmp(Indices{j},'lag') && MaxLag)||(strcmp(Indices{j},'lead') && MaxLead)
        [junk,row,col]=find(lead_lag_incidence(iter,:));
        eval(['A',Indices{j},'(:,row)=jacobia_(:,col);'])
        iter=iter+1;
    end
end
B=jacobia_(:,nnz(iyv)+1:end);

%%% MAIN ENGINE %%%
qz_criterium = options_.qz_criterium;
solve_maxit = options_.dp.maxit;
discretion_tol = options_.discretionary_tol;

if ~isempty(Hold)
    [H,G,info]=discretionary_policy_engine(Alag,A0,Alead,B,W,instr_id,beta,solve_maxit,discretion_tol,qz_criterium,Hold);
else
    [H,G,info]=discretionary_policy_engine(Alag,A0,Alead,B,W,instr_id,beta,solve_maxit,discretion_tol,qz_criterium);
end

if info
    dr=[];
    return
else
    Hold=H; %save previous solution
            % Hold=[]; use this line if persistent command is not used.
end
% set the state
dr=oo_.dr;
dr.ys =zeros(endo_nbr,1);
dr=set_state_space(dr,M_,options_);
order_var=dr.order_var;

T=H(order_var,order_var);
dr.ghu=G(order_var,:);
Selection=lead_lag_incidence(1,order_var)>0;%select state variables
dr.ghx=T(:,Selection);

ys=NondistortionarySteadyState(M_);
dr.ys=ys; % <--- dr.ys =zeros(NewEndo_nbr,1);

function ys=NondistortionarySteadyState(M_)
if exist([M_.fname,'_steadystate.m'],'file')
    eval(['ys=',M_.fname,'_steadystate.m;'])
else
    ys=zeros(M_.endo_nbr,1);
end
