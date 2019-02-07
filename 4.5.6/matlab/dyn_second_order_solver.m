function dr = dyn_second_order_solver(jacobia,hessian_mat,dr,M_,threads_ABC,threads_BC)

%@info:
%! @deftypefn {Function File} {@var{dr} =} dyn_second_order_solver (@var{jacobia},@var{hessian_mat},@var{dr},@var{M_},@var{threads_ABC},@var{threads_BC})
%! @anchor{dyn_second_order_solver}
%! @sp 1
%! Computes the second order reduced form of the DSGE model
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item jacobia
%! Matrix containing the Jacobian of the model
%! @item hessian_mat
%! Matrix containing the second order derivatives of the model
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @item M_
%! Matlab's structure describing the model (initialized by @code{dynare}).
%! @item threads_ABC
%! Integer controlling number of threads in A_times_B_kronecker_C
%! @item threads_BC
%! Integer controlling number of threads in sparse_hessian_times_B_kronecker_C
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @end table
%! @end deftypefn
%@eod:

% Copyright (C) 2001-2017 Dynare Team
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

dr.ghxx = [];
dr.ghuu = [];
dr.ghxu = [];
dr.ghs2 = [];
Gy = dr.Gy;

kstate = dr.kstate;
nstatic = M_.nstatic;
nfwrd = M_.nfwrd;
nspred = M_.nspred;
nboth = M_.nboth;
nsfwrd = M_.nsfwrd;
order_var = dr.order_var;
nd = size(kstate,1);
lead_lag_incidence = M_.lead_lag_incidence;

np = nd - nsfwrd;

k1 = nonzeros(lead_lag_incidence(:,order_var)');
kk = [k1; length(k1)+(1:M_.exo_nbr+M_.exo_det_nbr)'];
nk = size(kk,1);
kk1 = reshape([1:nk^2],nk,nk);
kk1 = kk1(kk,kk);
% reordering second order derivatives
hessian_mat = hessian_mat(:,kk1(:));

zx = zeros(np,np);
zu=zeros(np,M_.exo_nbr);
zx(1:np,:)=eye(np);
k0 = [1:M_.endo_nbr];
gx1 = dr.ghx;
hu = dr.ghu(nstatic+[1:nspred],:);
k0 = find(lead_lag_incidence(M_.maximum_endo_lag+1,order_var)');
zx = [zx; gx1(k0,:)];
zu = [zu; dr.ghu(k0,:)];
k1 = find(lead_lag_incidence(M_.maximum_endo_lag+2,order_var)');
zu = [zu; gx1(k1,:)*hu];
zx = [zx; gx1(k1,:)*Gy];
zx=[zx; zeros(M_.exo_nbr,np);zeros(M_.exo_det_nbr,np)];
zu=[zu; eye(M_.exo_nbr);zeros(M_.exo_det_nbr,M_.exo_nbr)];
[nrzx,nczx] = size(zx);

[rhs, err] = sparse_hessian_times_B_kronecker_C(hessian_mat,zx,threads_BC);
mexErrCheck('sparse_hessian_times_B_kronecker_C', err);
rhs = -rhs;

%lhs
n = M_.endo_nbr+sum(kstate(:,2) > M_.maximum_endo_lag+1 & kstate(:,2) < M_.maximum_endo_lag+M_.maximum_endo_lead+1);
A = zeros(M_.endo_nbr,M_.endo_nbr);
B = zeros(M_.endo_nbr,M_.endo_nbr);
A(:,k0) = jacobia(:,nonzeros(lead_lag_incidence(M_.maximum_endo_lag+1,order_var)));
% variables with the highest lead
k1 = find(kstate(:,2) == M_.maximum_endo_lag+2);
% Jacobian with respect to the variables with the highest lead
fyp = jacobia(:,kstate(k1,3)+nnz(M_.lead_lag_incidence(M_.maximum_endo_lag+1,:)));
B(:,nstatic+M_.npred+1:end) = fyp;
[junk,k1,k2] = find(M_.lead_lag_incidence(M_.maximum_endo_lag+M_.maximum_endo_lead+1,order_var));
A(1:M_.endo_nbr,nstatic+1:nstatic+nspred)=...
    A(1:M_.endo_nbr,nstatic+[1:nspred])+fyp*gx1(k1,1:nspred);
C = Gy;
D = [rhs; zeros(n-M_.endo_nbr,size(rhs,2))];


[err, dr.ghxx] = gensylv(2,A,B,C,D);
mexErrCheck('gensylv', err);

%ghxu
%rhs
hu = dr.ghu(nstatic+1:nstatic+nspred,:);
[rhs, err] = sparse_hessian_times_B_kronecker_C(hessian_mat,zx,zu,threads_BC);
mexErrCheck('sparse_hessian_times_B_kronecker_C', err);

hu1 = [hu;zeros(np-nspred,M_.exo_nbr)];
[nrhx,nchx] = size(Gy);
[nrhu1,nchu1] = size(hu1);

[abcOut,err] = A_times_B_kronecker_C(dr.ghxx,Gy,hu1,threads_ABC);
mexErrCheck('A_times_B_kronecker_C', err);
B1 = B*abcOut;
rhs = -[rhs; zeros(n-M_.endo_nbr,size(rhs,2))]-B1;


%lhs
dr.ghxu = A\rhs;

%ghuu
%rhs
[rhs, err] = sparse_hessian_times_B_kronecker_C(hessian_mat,zu,threads_BC);
mexErrCheck('sparse_hessian_times_B_kronecker_C', err);

[B1, err] = A_times_B_kronecker_C(B*dr.ghxx,hu1,threads_ABC);
mexErrCheck('A_times_B_kronecker_C', err);
rhs = -[rhs; zeros(n-M_.endo_nbr,size(rhs,2))]-B1;

%lhs
dr.ghuu = A\rhs;

% dr.ghs2
% derivatives of F with respect to forward variables
% reordering predetermined variables in diminishing lag order
O1 = zeros(M_.endo_nbr,nstatic);
O2 = zeros(M_.endo_nbr,M_.endo_nbr-nstatic-nspred);
LHS = zeros(M_.endo_nbr,M_.endo_nbr);
LHS(:,k0) = jacobia(:,nonzeros(lead_lag_incidence(M_.maximum_endo_lag+1,order_var)));
RHS = zeros(M_.endo_nbr,M_.exo_nbr^2);
gu = dr.ghu;
guu = dr.ghuu;
E = eye(M_.endo_nbr);
kh = reshape([1:nk^2],nk,nk);
kp = sum(kstate(:,2) <= M_.maximum_endo_lag+1);
E1 = [eye(nspred); zeros(kp-nspred,nspred)];
H = E1;
hxx = dr.ghxx(nstatic+[1:nspred],:);
[junk,k2a,k2] = find(M_.lead_lag_incidence(M_.maximum_endo_lag+2,order_var));
k3 = nnz(M_.lead_lag_incidence(1:M_.maximum_endo_lag+1,:))+(1:M_.nsfwrd)';
[B1, err] = sparse_hessian_times_B_kronecker_C(hessian_mat(:,kh(k3,k3)),gu(k2a,:),threads_BC);
mexErrCheck('sparse_hessian_times_B_kronecker_C', err);
RHS = RHS + jacobia(:,k2)*guu(k2a,:)+B1;

% LHS
LHS = LHS + jacobia(:,k2)*(E(k2a,:)+[O1(k2a,:) dr.ghx(k2a,:)*H O2(k2a,:)]);

RHS = RHS*M_.Sigma_e(:);
dr.fuu = RHS;
%RHS = -RHS-dr.fbias;
RHS = -RHS;
dr.ghs2 = LHS\RHS;

% deterministic exogenous variables
if M_.exo_det_nbr > 0
end
