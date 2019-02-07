function [dr,info] = dyn_risky_steadystate_solver(ys0,M, ...
                                                  dr,options,oo)

%@info:
%! @deftypefn {Function File} {[@var{dr},@var{info}] =} dyn_risky_steadystate_solver (@var{ys0},@var{M},@var{dr},@var{options},@var{oo})
%! @anchor{dyn_risky_steadystate_solver}
%! @sp 1
%! Computes the second order risky steady state and first and second order reduced form of the DSGE model.
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item ys0
%! Vector containing a guess value for the risky steady state
%! @item M
%! Matlab's structure describing the model (initialized by @code{dynare}).
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @item options
%! Matlab's structure describing the options (initialized by @code{dynare}).
%! @item oo
%! Matlab's structure gathering the results (initialized by @code{dynare}).
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @ @var
%! @item dr
%! Matlab's structure describing the reduced form solution of the model.
%! @item info
%! Integer scalar, error code.
%! @sp 1
%! @table @ @code
%! @item info==0
%! No error.
%! @item info==1
%! The model doesn't determine the current variables uniquely.
%! @item info==2
%! MJDGGES returned an error code.
%! @item info==3
%! Blanchard & Kahn conditions are not satisfied: no stable equilibrium.
%! @item info==4
%! Blanchard & Kahn conditions are not satisfied: indeterminacy.
%! @item info==5
%! Blanchard & Kahn conditions are not satisfied: indeterminacy due to rank failure.
%! @item info==6
%! The jacobian evaluated at the deterministic steady state is complex.
%! @item info==19
%! The steadystate routine thrown an exception (inconsistent deep parameters).
%! @item info==20
%! Cannot find the steady state, info(2) contains the sum of square residuals (of the static equations).
%! @item info==21
%! The steady state is complex, info(2) contains the sum of square of imaginary parts of the steady state.
%! @item info==22
%! The steady has NaNs.
%! @item info==23
%! M_.params has been updated in the steadystate routine and has complex valued scalars.
%! @item info==24
%! M_.params has been updated in the steadystate routine and has some NaNs.
%! @end table
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


info = 0;
lead_lag_incidence = M.lead_lag_incidence;
order_var = dr.order_var;
endo_nbr = M.endo_nbr;
exo_nbr = M.exo_nbr;

M.var_order_endo_names = M.endo_names(dr.order_var,:);

[junk,dr.i_fwrd_g,i_fwrd_f] = find(lead_lag_incidence(3,order_var));
dr.i_fwrd_f = i_fwrd_f;
nd = nnz(lead_lag_incidence) + M.exo_nbr;
dr.nd = nd;
kk = reshape(1:nd^2,nd,nd);
kkk = reshape(1:nd^3,nd^2,nd);
dr.i_fwrd2_f = kk(i_fwrd_f,i_fwrd_f);
dr.i_fwrd2a_f = kk(i_fwrd_f,:);
dr.i_fwrd3_f = kkk(dr.i_fwrd2_f,:);
dr.i_uu = kk(end-exo_nbr+1:end,end-exo_nbr+1:end);
if options.k_order_solver
    func = @risky_residuals_k_order;
else
    func = @risky_residuals;
end

if isfield(options,'portfolio') && options.portfolio == 1
    pm = portfolio_model_structure(M,options);

    x0 = ys0(pm.v_p);
    n = length(x0);
    [x, info] = solve1(@risky_residuals_ds,x0,1:n,1:n,0,options.gstep, ...
                       options.solve_tolf,options.solve_tolx, ...
                       options.steady.maxit,options.debug,pm,M,dr, ...
                       options,oo);
    if info
        error('DS approach can''t be computed')
    end
    %[x, info] = csolve(@risky_residuals_ds,x0,[],1e-10,100,M,dr,options,oo);
    %        ys0(l_var) = x;
    [resids,dr1] = risky_residuals_ds(x,pm,M,dr,options,oo);
    ys1 = dr1.ys;
else
    pm = model_structure(M,options);
end

[ys, info] = solve1(func,ys0,1:endo_nbr,1:endo_nbr,0,options.gstep, ...
                    options.solve_tolf,options.solve_tolx, ...
                    options.steady.maxit,options.debug,pm,M,dr,options,oo);
%    [ys, info] = csolve(func,ys0,[],1e-10,100,M,dr,options,oo);
if info
    error('RSS approach can''t be computed')
end
dr.ys = ys;

[resid,dr] = func(ys,pm,M,dr,options,oo);
dr.ghs2 = zeros(M.endo_nbr,1);

for i=1:M.endo_nbr
    if isfield(options,'portfolio') && options.portfolio == 1
        disp(sprintf('%16s %12.6f %12.6f',M.endo_names(i,:),ys1(i), ...
                     ys(i)))
    else
        disp(sprintf('%16s %12.6f %12.6f',M.endo_names(i,:),ys(i)))
    end
end

end

function [resid,dr] = risky_residuals(ys,pm,M,dr,options,oo)

lead_lag_incidence = M.lead_lag_incidence;
iyv = lead_lag_incidence';
iyv = iyv(:);
iyr0 = find(iyv) ;

if M.exo_nbr == 0
    oo.exo_steady_state = [] ;
end

z = repmat(ys,1,3);
z = z(iyr0) ;
[resid1,d1,d2] = feval([M.fname '_dynamic'],z,...
                       [oo.exo_simul ...
                    oo.exo_det_simul], M.params, dr.ys, 2);
if ~isreal(d1) || ~isreal(d2)
    pause
end

if options.use_dll
    % In USE_DLL mode, the hessian is in the 3-column sparse representation
    d2 = sparse(d2(:,1), d2(:,2), d2(:,3), ...
                size(d1, 1), size(d1, 2)*size(d1, 2));
end

if isfield(options,'portfolio') && options.portfolio == 1
    pm = portfolio_model_structure(M,options);
    x = ys(pm.v_p);
    dr = first_step_ds(x,pm,M,dr,options,oo);
    dr.ys = ys;
else
    pm = model_structure(M,options);
    [dr,info] = dyn_first_order_solver(d1,M,dr,options,0);
    if info
        print_info(info,options.noprint,options);
    end
    dr = dyn_second_order_solver(d1,d2,dr,M,...
                                 options.threads.kronecker.A_times_B_kronecker_C,...
                                 options.threads.kronecker.sparse_hessian_times_B_kronecker_C);
end

gu1 = dr.ghu(pm.i_fwrd_g,:);

resid = resid1+0.5*(d1(:,pm.i_fwrd_f1)*dr.ghuu(pm.i_fwrd_g,:)+ ...
                    d2(:,pm.i_fwrd_f2)*kron(gu1,gu1))*vec(M.Sigma_e);
end

function [resid,dr] = risky_residuals_ds(x,pm,M,dr,options,oo)

v_p = pm.v_p;
v_np = pm.v_np;

% computing steady state of non-portfolio variables  consistent with
% assumed portfolio
dr.ys(v_p) = x;
ys0 = dr.ys(v_np);
f_h =str2func([M.fname '_static']);
[dr.ys(v_np),info] = csolve(@ds_static_model,ys0,[],1e-10,100,f_h,x,pm.eq_np,v_np,v_p, ...
                            M.endo_nbr,M.exo_nbr,M.params);
if info
    error('can''t compute non-portfolio steady state')
end

dr_np = first_step_ds(x,pm,M,dr,options,oo);

lead_lag_incidence = M.lead_lag_incidence;
iyv = lead_lag_incidence';
iyv = iyv(:);
iyr0 = find(iyv) ;

z = repmat(dr.ys,1,3);
z = z(iyr0) ;
[resid1,d1,d2] = feval([M.fname '_dynamic'],z,...
                       [oo.exo_simul ...
                    oo.exo_det_simul], M.params, dr.ys, 2);
if ~isreal(d1) || ~isreal(d2)
    pause
end

if options.use_dll
    % In USE_DLL mode, the hessian is in the 3-column sparse representation
    d2 = sparse(d2(:,1), d2(:,2), d2(:,3), ...
                size(d1, 1), size(d1, 2)*size(d1, 2));
end


gu1 = dr_np.ghu(pm.i_fwrd_g,:);

resid = resid1+0.5*(d2(:,pm.i_fwrd_f2)*kron(gu1,gu1))*vec(M.Sigma_e);

resid = resid(pm.eq_p)
end

function dr_np = first_step_ds(x,pm,M,dr,options,oo)

lead_lag_incidence = M.lead_lag_incidence;
iyv = lead_lag_incidence';
iyv = iyv(:);
iyr0 = find(iyv) ;

ys = dr.ys;
ys(pm.v_p) = x;

z = repmat(ys,1,3);
z = z(iyr0) ;
[resid1,d1,d2] = feval([M.fname '_dynamic'],z,...
                       [oo.exo_simul ...
                    oo.exo_det_simul], M.params, dr.ys, 2);
if ~isreal(d1) || ~isreal(d2)
    pause
end

if options.use_dll
    % In USE_DLL mode, the hessian is in the 3-column sparse representation
    d2 = sparse(d2(:,1), d2(:,2), d2(:,3), ...
                size(d1, 1), size(d1, 2)*size(d1, 2));
end

d1_np = d1(pm.eq_np,pm.i_d1_np);
d2_np = d2(pm.eq_np,pm.i_d2_np);

[dr_np,info] = dyn_first_order_solver(d1_np,pm.M_np,pm.dr_np,options,0);
if info
    print_info(info, 0, options);
    return
end

dr_np = dyn_second_order_solver(d1_np,d2_np,dr_np,pm.M_np,...
                                options.threads.kronecker.A_times_B_kronecker_C,...
                                options.threads.kronecker.sparse_hessian_times_B_kronecker_C);
end

function [resid,dr] = risky_residuals_k_order(ys,pm,M,dr,options,oo)
exo_nbr = M.exo_nbr;
endo_nbr = M.endo_nbr;

iyv = M.lead_lag_incidence';
iyv = iyv(:);
iyr0 = find(iyv) ;

if exo_nbr == 0
    oo.exo_steady_state = [] ;
end

z = repmat(ys,1,3);
z = z(iyr0) ;
[resid1,d1,d2] = feval([M.fname '_dynamic'],z,...
                       [oo.exo_simul ...
                    oo.exo_det_simul], M.params, dr.ys, 2);

if isfield(options,'portfolio') && options.portfolio == 1
    eq_np = pm.eq_np;

    d1_np = d1(eq_np,pm.i_d1_np);
    d2_np = d2(eq_np,pm.i_d2_np);

    M_np = pm.M_np;
    dr_np = pm.dr_np;

    [dr_np,info] = dyn_first_order_solver(d1_np,pm.M_np,pm.dr_np,options,0);
    if info
        print_info(info, 0, options);
        return
    end

    dr_np = dyn_second_order_solver(d1_np,d2_np,dr_np,pm.M_np,...
                                    options.threads.kronecker.A_times_B_kronecker_C,...
                                    options.threads.kronecker.sparse_hessian_times_B_kronecker_C);
end

i_fwrd_f1 = pm.i_fwrd_f1;
i_fwrd_f2 = pm.i_fwrd_f2;
i_fwrd_f3 = pm.i_fwrd_f3;
i_fwrd_g = pm.i_fwrd_g;
gu1 = dr_np.ghu(i_fwrd_g,:);
ghuu = dr_np.ghuu;

resid = resid1+0.5*(d1(:,i_fwrd_f1)*ghuu(i_fwrd_g,:)+d2(:,i_fwrd_f2)* ...
                    kron(gu1,gu1))*vec(M.Sigma_e);

if nargout > 1
    [resid1,d1,d2,d3] = feval([M.fname '_dynamic'],z,...
                              [oo.exo_simul ...
                        oo.exo_det_simul], M.params, dr.ys, 2);


    [a,b,c] = find(d2(eq_np,pm.i_d2_np));
    d2_np = [a b c];

    [a,b,c] = find(d3(eq_np,pm.i_d3_np));
    d3_np = [a b c];

    options.order = 3;
    % space holder, unused by k_order_pertrubation
    dr_np.ys = dr.ys(pm.v_np);
    nu2 = exo_nbr*(exo_nbr+1)/2;
    nu3 = exo_nbr*(exo_nbr+1)*(exo_nbr+2)/3;
    M_np.NZZDerivatives = [nnz(d1_np); nnz(d2_np); nnz(d3_np)];
    [err,g_0, g_1, g_2, g_3] = k_order_perturbation(dr_np,M_np,options,d1_np,d2_np,d3_np);
    mexErrCheck('k_order_perturbation', err);

    gu1 = g_1(i_fwrd_g,end-exo_nbr+1:end);
    ghuu = unfold2(g_2(:,end-nu2+1:end),exo_nbr);
    ghsuu = get_ghsuu(g_3,size(g_1,2),exo_nbr);

    i_fwrd1_f2 = pm.i_fwrd1_f2;
    i_fwrd1_f3 = pm.i_fwrd1_f3;
    n = size(d1,2);
    d1b = d1 + 0.5*( ...
        d1(:,i_fwrd_f1)*...
        d2(:,i_fwrd1_f2)*kron(eye(n),dr_np.ghuu(i_fwrd_g,:)*vec(M.Sigma_e))...
        + 0.5*d3(:,i_fwrd1_f3)*kron(eye(n),kron(gu1,gu1)*vec(M.Sigma_e)));
    format short
    kk1 = [nonzeros(M.lead_lag_incidence(:,1:6)'); ...
           nnz(M.lead_lag_incidence)+[1; 2]]
    kk2 = [nonzeros(M.lead_lag_incidence(:,1:6)'); ...
           nnz(M.lead_lag_incidence)+[3; 4]]
    format short
    gu1
    kron(gu1,gu1)*vec(M.Sigma_e)
    disp(d1(:,:))
    disp(d1b(:,:))
    aa2=d2(:,i_fwrd1_f2)*kron(eye(n),dr_np.ghuu(i_fwrd_g,:)*vec(M.Sigma_e));
    aa3=d3(:,i_fwrd1_f3)*kron(eye(n),kron(gu1,gu1)*vec(M.Sigma_e));
    disp(d3(4,7+6*n+6*n*n))
    disp(d3(4,8+16*n+17*n*n))   %8,17,18
    disp(d3(4,8+17*n+16*n*n))   %8,17,18
    disp(d3(4,7*n+17+17*n*n))   %8,17,18
    disp(d3(4,7*n+18+16*n*n))   %8,17,18
    disp(d3(4,7*n*n+16*n+18))   %8,17,18
    disp(d3(4,7*n*n+17+17*n))   %8,17,18
    pause
    disp(aa2(:,kk1))
    disp(aa2(:,kk2))
    disp(aa3(:,kk1))
    disp(aa3(:,kk2))
    [dr,info] = dyn_first_order_solver(d1b,M,dr,options,0);
    if info
        print_info(info, 0, options);
        return
    end

    disp_dr(dr,dr.order_var,[]);

end
end

function y=unfold2(x,n)
y = zeros(size(x,1),n*n);
k = 1;
for i=1:n
    for j=i:n
        y(:,(i-1)*n+j) = x(:,k);
        if i ~= j
            y(:,(j-1)*n+i) = x(:,k);
        end
        k = k+1;
    end
end
end

function y=unfold3(x,n)
y = zeros(size(x,1),n*n*n);
k = 1;
for i=1:n
    for j=i:n
        for m=j:n
            y(:,(i-1)*n*n+(j-1)*n+m) = x(:,k);
            y(:,(i-1)*n*n+(m-1)*n+j) = x(:,k);
            y(:,(j-1)*n*n+(i-1)*n+m) = x(:,k);
            y(:,(j-1)*n*n+(m-1)*n+i) = x(:,k);
            y(:,(m-1)*n*n+(i-1)*n+j) = x(:,k);
            y(:,(m-1)*n*n+(j-1)*n+i) = x(:,k);

            k = k+1;
        end
    end
end
end

function pm  = model_structure(M,options)


lead_index = M.maximum_endo_lag+2;
lead_lag_incidence = M.lead_lag_incidence;
dr = struct();
dr = set_state_space(dr,M,options);
pm.i_fwrd_g = find(lead_lag_incidence(lead_index,dr.order_var)');

i_fwrd_f1 = nonzeros(lead_lag_incidence(lead_index,dr.order_var));
pm.i_fwrd_f1 = i_fwrd_f1;
n = nnz(lead_lag_incidence)+M.exo_nbr;
ih = reshape(1:n*n,n,n);
i_fwrd_f2 = ih(i_fwrd_f1,i_fwrd_f1);
pm.i_fwrd_f2 = i_fwrd_f2(:);
i_fwrd1_f2 = ih(i_fwrd_f1,:);
pm.i_fwrd1_f2 = i_fwrd1_f2(:);

ih = reshape(1:n*n*n,n,n,n);
i_fwrd_f3 = ih(i_fwrd_f1,i_fwrd_f1,i_fwrd_f1);
pm.i_fwrd_f3 = i_fwrd_f3(:);
i_fwrd1_f3 = ih(i_fwrd_f1,i_fwrd_f1,:);
pm.i_fwrd1_f3 = i_fwrd1_f3(:);
end

function pm  = portfolio_model_structure(M,options)

i_d3_np = [];
i_d3_p = [];

lead_index = M.maximum_endo_lag+2;
lead_lag_incidence = M.lead_lag_incidence;
eq_tags = M.equations_tags;
n_tags = size(eq_tags,1);
eq_p = cell2mat(eq_tags(strcmp(eq_tags(:,2), ...
                               'portfolio'),1));
pm.eq_p = eq_p;
pm.eq_np = setdiff(1:M.endo_nbr,eq_p);
v_p = zeros(n_tags,1);
for i=1:n_tags
    v_p(i) = find(strncmp(eq_tags(i,3),M.endo_names, ...
                          length(cell2mat(eq_tags(i,3)))));
end
if any(lead_lag_incidence(lead_index,v_p))
    error(['portfolio variables appear in the model as forward ' ...
           'variable'])
end
pm.v_p = v_p;
v_np = setdiff(1:M.endo_nbr,v_p);
pm.v_np = v_np;
lli_np = lead_lag_incidence(:,v_np)';
k = find(lli_np);
lead_lag_incidence_np = lli_np;
lead_lag_incidence_np(k) = 1:nnz(lli_np);
lead_lag_incidence_np = lead_lag_incidence_np';
pm.lead_lag_incidence_np = lead_lag_incidence_np;
i_d1_np = [nonzeros(lli_np); nnz(lead_lag_incidence)+(1:M.exo_nbr)'];
pm.i_d1_np = i_d1_np;

n = nnz(lead_lag_incidence)+M.exo_nbr;
ih = reshape(1:n*n,n,n);
i_d2_np = ih(i_d1_np,i_d1_np);
pm.i_d2_np = i_d2_np(:);

ih = reshape(1:n*n*n,n,n,n);
i_d3_np = ih(i_d1_np,i_d1_np,i_d1_np);
pm.i_d3_np = i_d3_np(:);

M_np = M;
M_np.lead_lag_incidence = lead_lag_incidence_np;
M_np.lead_lag_incidence = lead_lag_incidence_np;
M_np.endo_nbr = length(v_np);
M_np.endo_names = M.endo_names(v_np,:);
dr_np = struct();
dr_np = set_state_space(dr_np,M_np,options);
pm.dr_np = dr_np;
M_np.var_order_endo_names = M_np.endo_names(dr_np.order_var,:);
pm.M_np = M_np;
pm.i_fwrd_g = find(lead_lag_incidence_np(lead_index,dr_np.order_var)');

i_fwrd_f1 = nonzeros(lead_lag_incidence(lead_index,:));
pm.i_fwrd_f1 = i_fwrd_f1;
n = nnz(lead_lag_incidence)+M.exo_nbr;
ih = reshape(1:n*n,n,n);
i_fwrd_f2 = ih(i_fwrd_f1,i_fwrd_f1);
pm.i_fwrd_f2 = i_fwrd_f2(:);
i_fwrd1_f2 = ih(i_fwrd_f1,:);
pm.i_fwrd1_f2 = i_fwrd1_f2(:);

ih = reshape(1:n*n*n,n,n,n);
i_fwrd_f3 = ih(i_fwrd_f1,i_fwrd_f1,i_fwrd_f1);
pm.i_fwrd_f3 = i_fwrd_f3(:);
i_fwrd1_f3 = ih(i_fwrd_f1,i_fwrd_f1,:);
pm.i_fwrd1_f3 = i_fwrd1_f3(:);
end

function r=ds_static_model(y0,f_h,p0,eq_np,v_np,v_p,endo_nbr,exo_nbr,params)
ys = zeros(endo_nbr,1);
ys(v_p) = p0;
ys(v_np) = y0;
r = f_h(ys,zeros(exo_nbr,1),params);
r = r(eq_np);
end

function ghsuu = get_ghsuu(g,ns,nx)
nxx = nx*(nx+1)/2;
m1 = 0;
m2 = ns*(ns+1)/2;
kk = 1:(nx*nx);
ghsuu = zeros(size(g,1),(ns*nx*nx));

for i=1:n
    j = m1+(1:m2);
    k = j(end-nxx+1:end);
    ghsuu(:,kk) = unfold2(g(:,k),nx);
    m1 = m1+m2;
    m2 = m2 - (n-i+1);
    kk = kk + nx*nx;
end
end