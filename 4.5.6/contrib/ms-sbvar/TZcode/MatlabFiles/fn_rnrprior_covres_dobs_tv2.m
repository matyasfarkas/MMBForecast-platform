function [Pi_bar,H0tldcell_inv,Hptldcell_inv] ...
                      = fn_rnrprior_covres_dobs_tv2(nvar,nStates,indxScaleStates,q_m,lags,xdgel,mu,indxDummy,Ui,Vi,hpmsmd,indxmsmdeqn,nexo,asym0,asymp)
% Differs from fn_rnrprior_covres_dobs(): linear restrictions (Ui and Vi) have been incorported in fn_rnrprior_covres_dobs_tv?().
% Differs from fn_rnrprior_covres_dobs_tv():  allows an option to scale up the prior variance by nStates or not scale at all,
%                so that the prior value is the same as the constant VAR when the parameters in all states are the same.
%
% Only works for the nexo=1 (constant term) case.  To extend this to other exogenous variables, see fn_dataxy.m.  01/14/03.
% Differs from fn_rnrprior_covres_tv.m in that dummy observations are included as part of the explicit prior.  See Forcast II, pp.68-69b.
% Exports random Bayesian prior of Sims and Zha with asymmetric rior with linear restrictions already applied
%   and with dummy observations (i.e., mu(5) and mu(6)) used as part of an explicit prior.
% This function allows for prior covariances for the MS and MD equations to achieve liquidity effects.
% See Waggoner and Zha's Gibbs sampling paper and TVBVAR NOTES pp. 71k.0 and 50-61.
%
% nvar:  number of endogenous variables
% nStates:  Number of states.
% indxScaleStates: if 0, no scale adjustment in the prior variance for the number of states in the function fn_rnrprior_covres_dobs_tv2();
%                  if 1: allows a scale adjustment, marking the prior variance bigger by the number of states.
% q_m:  quarter or month
% lags: the maximum length of lag
% xdgel: T*nvar endogenous-variable matrix of raw or original data (no manipulation involved) with sample size including lags.
%       Order of columns: (1) nvar endogenous variables; (2) constants will be automatically put in the last column.
%       Used only to get variances of residuals for mu(1)-mu(5) and for dummy observations mu(5) and mu(6).
% mu: 6-by-1 vector of hyperparameters (the following numbers for Atlanta Fed's forecast), where
%          mu(5) and mu(6) are NOT used here.  See fn_dataxy.m for using mu(5) and mu(6).
%       mu(1): overall tightness and also for A0;  (0.57)
%       mu(2): relative tightness for A+;  (0.13)
%       mu(3): relative tightness for the constant term;  (0.1).  NOTE: for other
%               exogenous terms, the variance of each exogenous term must be taken into
%               acount to eliminate the scaling factor.
%       mu(4): tightness on lag decay;  (1)
%       mu(5): weight on nvar sums of coeffs dummy observations (unit roots);  (5)
%       mu(6): weight on single dummy initial observation including constant
%               (cointegration, unit roots, and stationarity);  (5)
%       NOTE: for this function, mu(5) and mu(6) are not used.  See fn_dataxy.m for using mu(5) and mu(6).
% indxDummy:  1: uses dummy observations to form part of an explicit prior; 0: no dummy observations as part of the prior.
% Ui: nvar-by-1 cell.  In each cell, nvar-by-qi*si orthonormal basis for the null of the ith
%           equation contemporaneous restriction matrix where qi is the number of free parameters
%           within the state and si is the number of free states.
%           With this transformation, we have ai = Ui*bi or Ui'*ai = bi where ai is a vector
%           of total original parameters and bi is a vector of free parameters. When no
%           restrictions are imposed, we have Ui = I.  There must be at least one free
%           parameter left for the ith equation in the order of [a_i for 1st state, ..., a_i for last state].
% Vi: nvar-by-1 cell.  In each cell, k-by-ri*ti orthonormal basis for the null of the ith
%           equation lagged restriction matrix where k is a total of exogenous variables and
%           ri is the number of free parameters within the state and ti is the number of free states.
%           With this transformation, we have fi = Vi*gi
%           or Vi'*fi = gi where fi is a vector of total original parameters and gi is a
%           vector of free parameters.  The ith equation is in the order of [nvar variables
%           for 1st lag and 1st state, ..., nvar variables for last lag and 1st state, const for 1st state, nvar
%           variables for 1st lag and 2nd state, nvar variables for last lag and 2nd state, const for 2nd state, and so on].
% hpmsmd: 2-by-1 hyperparameters with -1<h1=hpmsmd(1)<=0 for the MS equation and 0<=h2=hpmsmd(2)<1 the MD equation.  Consider a1*R + a2*M.
%          The term h1*var(a1)*var(a2) is the prior covariance of a1 and a2 for MS, equivalent to penalizing the same sign of a1 and a2.
%          The term h2*var(a1)*var(a2) is the prior covariance of a1 and a2 for MD, equivalent to penalizing opposite signs of a1 and a2.
%          This will give us a liquidity effect.  If hpmsmd=0, no such restrictions will be imposed.
% indxmsmdeqn: 4-by-1 index for the locations of the MS and MD equation and for the locations of M and R.
%               indxmsmdeqn(1) for MS and indxmsmdeqn(2) for MD.
%               indxmsmdeqn(3) for M and indxmsmdeqn(4) for R.
% nexo:  number of exogenous variables (if not specified, nexo=1 (constant) by default).
%         The constant term is always put to the last of all endogenous and exogenous variables.
% asym0: nvar-by-nvar asymmetric prior on A0.  Column -- equation.
%        If ones(nvar,nvar), symmetric prior;  if not, relative (asymmetric) tightness on A0.
% asymp: ncoef-1-by-nvar asymmetric prior on A+ bar constant.  Column -- equation.
%        If ones(ncoef-1,nvar), symmetric prior;  if not, relative (asymmetric) tightness on A+.
% --------------------
% Pi_bar: ncoef-by-nvar matrix for the ith equation under random walk.  Same for all equations
% H0tldcell_inv: cell(nvar,1).  The ith cell represents the ith equation, where the dim is
%         qi*si-by-qi*si.  The inverse of H0tld on p.60.
% Hptldcell_inv: cell(nvar,1).  The ith cell represents the ith equation, where the dim is
%         ri*ti-by-ri*ti.The inverse of Hptld on p.60.
%
% Differs from fn_rnrprior_covres_dobs(): linear restrictions (Ui and Vi) have been incorported in fn_rnrprior_covres_dobs_tv?().
% Differs from fn_rnrprior_covres_dobs_tv():  allows an option to scale up the prior variance by nStates or not scale at all.
%                so that the prior value is the same as the constant VAR when the parameters in all states are the same.
% Tao Zha, February 2000.  Revised, September 2000, 2001, February, May 2003, May 2004.
%
% Copyright (C) 1997-2012 Tao Zha
%
% This free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% It is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% If you did not received a copy of the GNU General Public License
% with this software, see <http://www.gnu.org/licenses/>.
%


if nargin<=12, nexo=1; end   % <<>>1
ncoef = nvar*lags+nexo;  % Number of coefficients in *each* equation for each state, RHS coefficients only.
ncoefsts = nStates*ncoef;  % Number of coefficients in *each* equation in all states, RHS coefficients only.

H0tldcell_inv=cell(nvar,1);  % inv(H0tilde) for different equations under asymmetric prior.
Hptldcell_inv=cell(nvar,1);  % inv(H+tilde) for different equations under asymmetric prior.

%*** Constructing Pi_bar for the ith equation under the random walk assumption
Pi_bar = zeros(ncoef,nvar);   % same for all equations
Pi_bar(1:nvar,1:nvar) = eye(nvar);   % random walk

%
%@@@ Prepared for Bayesian prior
%
%
% ** monthly lag decay in order to match quarterly decay: a*exp(bl) where
% **  l is the monthly lag.  Suppose quarterly decay is 1/x where x=1,2,3,4.
% **  Let the decay of l1 (a*exp(b*l1)) match that of x1 (say, beginning: 1/1)
% **  and the decay of l2 (a*exp(b*l2)) match that of x2 (say, end: 1/5),
% **  we can solve for a and b which are
% **      b = (log_x1-log_x2)/(l1-l2), and a = x1*exp(-b*l1).
if q_m==12
   l1 = 1;   % 1st month == 1st quarter
   xx1 = 1;   % 1st quarter
   l2 = lags;   % last month
   xx2 = 1/((ceil(lags/3))^mu(4));   % last quarter
   %xx2 = 1/6;   % last quarter
   % 3rd quarter:  i.e., we intend to let decay of the 6th month match
   %    that of the 3rd quarter, so that the 6th month decays a little
   %    faster than the second quarter which is 1/2.
   if lags==1
      b = 0;
   else
      b = (log(xx1)-log(xx2))/(l1-l2);
   end
   a = xx1*exp(-b*l1);
end



%
% *** specify the prior for each equation separately, SZ method,
% ** get the residuals from univariate regressions.
%
sgh = zeros(nvar,1);        % square root
sgsh = sgh;              % square
nSample=size(xdgel,1);  % sample size-lags
yu = xdgel;
C = ones(nSample,1);
for k=1:nvar
   [Bk,ek,junk1,junk2,junk3,junk4] = sye([yu(:,k) C],lags);
   clear Bk junk1 junk2 junk3 junk4;
   sgsh(k) = ek'*ek/(nSample-lags);
   sgh(k) = sqrt(sgsh(k));
end
% ** prior variance for A0(:,1), same for all equations!!!
sg0bid = zeros(nvar,1);  % Sigma0_bar diagonal only for the ith equation
for j=1:nvar
   sg0bid(j) = 1/sgsh(j);    % sgsh = sigmai^2
end
% ** prior variance for lagged and exogeous variables, same for all equations
sgpbid = zeros(ncoef,1);     % Sigma_plus_bar, diagonal, for the ith equation
for i = 1:lags
   if (q_m==12)
      lagdecay = a*exp(b*i*mu(4));
   end
   %
   for j = 1:nvar
      if (q_m==12)
         % exponential decay to match quarterly decay
         sgpbid((i-1)*nvar+j) = lagdecay^2/sgsh(j);  % ith equation
      elseif (q_m==4)
         sgpbid((i-1)*nvar+j) = (1/i^mu(4))^2/sgsh(j);  % ith equation
      else
			error('Incompatibility with lags, check the possible errors!!!')
         %warning('Incompatibility with lags, check the possible errors!!!')
         %return
      end
   end
end
%
if indxDummy   % Dummy observations as part of the explicit prior.
   ndobs=nvar+1;         % Number of dummy observations: nvar unit roots and 1 cointegration prior.
   phibar = zeros(ndobs,ncoef);
   %* constant term
   const = ones(nvar+1,1);
   const(1:nvar) = 0.0;
   phibar(:,ncoef) = const;      % the first nvar periods: no or zero constant!

   xdgelint = mean(xdgel(1:lags,:),1); % mean of the first lags initial conditions
   %* Dummies
   for k=1:nvar
      for m=1:lags
         phibar(ndobs,nvar*(m-1)+k) = xdgelint(k);
         phibar(k,nvar*(m-1)+k) = xdgelint(k);
         % <<>> multiply hyperparameter later
      end
   end
   phibar(1:nvar,:) = 1*mu(5)*phibar(1:nvar,:);    % standard Sims and Zha prior
   phibar(ndobs,:) = mu(6)*phibar(ndobs,:);
   [phiq,phir]=qr(phibar,0);
   xtxbar=phir'*phir;   % phibar'*phibar.   ncoef-by-ncoef.  Reduced (not full) rank.  See Forcast II, pp.69-69b.
end





%=================================================
%   Computing the (prior) covariance matrix for the posterior of A0, no data yet
%=================================================
%
%
% ** set up the conditional prior variance sg0bi and sgpbi.
sg0bida = mu(1)^2*sg0bid;   % ith equation
sgpbida = mu(1)^2*mu(2)^2*sgpbid;
sgpbida(ncoef-nexo+1:ncoef) = mu(1)^2*mu(3)^2;
          %<<>> No scaling adjustment has been made for exogenous terms other than constant
sgppbd = sgpbida(nvar+1:ncoef);    % corresponding to A++, in a Sims-Zha paper

Hptd = zeros(ncoef);
Hptdi=Hptd;
Hptd(ncoef,ncoef)=sgppbd(ncoef-nvar);
Hptdinv(ncoef,ncoef)=1./sgppbd(ncoef-nvar);
             % condtional on A0i, H_plus_tilde


if nargin<14   % <<>>1 Default is no asymmetric information
   asym0 = ones(nvar,nvar);  % if not ones, then we have relative (asymmetric) tightness
   asymp = ones(ncoef-1,nvar);    % for A+.  Column -- equation
end

%**** Asymmetric Information
%asym0 = ones(nvar,nvar);  % if not ones, then we have relative (asymmetric) tightness
%asymp = ones(ncoef-1,nvar);    % pp: plus without constant.  Column -- equation
%>>>>>> B: asymmetric prior variance for asymp <<<<<<<<
%
%for i = 1:lags
%   rowif = (i-1)*nvar+1;
%   rowil = i*nvar;
%     idmatw0 = 0.5;   % weight assigned to idmat0 in the formation of asymp
%	if (i==1)
%     asymp(rowif:rowil,:)=(1-idmatw0)*ones(nvar)+idmatw0*idmat0;  % first lag
%		                 % note:  idmat1 is already transposed.  Column -- equation
%	else
%     %asymp(rowif:rowil,1:nvar) = (1-idmatw0)*ones(nvar)+idmatw0*idmat0;
%                % <<<<<<< toggle +
%                % Note: already transposed, since idmat0 is transposed.
%				     % Meaning: column implies equation
%     asymp(rowif:rowil,1:nvar) = ones(nvar);
%                % >>>>>>> toggle -
%	end
%end
%
%>>>>>> E: asymmetric prior variance for asymp <<<<<<<<


%=================================================
%   Computing the final covariance matrix (S1,...,Sm) for the prior of A0,
%      and final Gb=(G1,...,Gm) for A+ if asymmetric prior or for
%      B if symmetric prior for A+
%=================================================
%
for i = 1:nvar
   %------------------------------
   % Introduce prior information on which variables "belong" in various equations.
   % In this first trial, we just introduce this information here, in a model-specific way.
   % Eventually this info has to be passed parametricly.  In our first shot, we just damp down
   % all coefficients except those on the diagonal.

   %*** For A0
   factor0=asym0(:,i);

   sg0bd = sg0bida.*factor0;  %  Note, this only works for the prior variance Sg(i)
                      % of a0(i) being diagonal.  If the prior variance Sg(i) is not
                      % diagonal, we have to the inverse to get inv(Sg(i)).
   %sg0bdinv = 1./sg0bd;
   % *    unconditional variance on A0+
   H0td = diag(sg0bd);    % unconditional
   %=== Correlation in the MS equation to get a liquidity effect.
   if (i==indxmsmdeqn(1))
      H0td(indxmsmdeqn(3),indxmsmdeqn(4)) = hpmsmd(1)*sqrt(sg0bida(indxmsmdeqn(3))*sg0bida(indxmsmdeqn(4)));
      H0td(indxmsmdeqn(4),indxmsmdeqn(3)) = hpmsmd(1)*sqrt(sg0bida(indxmsmdeqn(3))*sg0bida(indxmsmdeqn(4)));
   elseif (i==indxmsmdeqn(2))
      H0td(indxmsmdeqn(3),indxmsmdeqn(4)) = hpmsmd(2)*sqrt(sg0bida(indxmsmdeqn(3))*sg0bida(indxmsmdeqn(4)));
      H0td(indxmsmdeqn(4),indxmsmdeqn(3)) = hpmsmd(2)*sqrt(sg0bida(indxmsmdeqn(3))*sg0bida(indxmsmdeqn(4)));
   end
   H0tdinv = inv(H0td);
   %H0tdinv = diag(sg0bdinv);
   %
   if indxScaleStates
      H0tldcell_inv{i}=(Ui{i}'*kron(eye(nStates),H0tdinv/nStates))*Ui{i};
   else
      H0tldcell_inv{i}=(Ui{i}'*kron(eye(nStates),H0tdinv))*Ui{i};
   end



   %*** For A+
   if ~(lags==0)  % For A1 to remain random walk properties
      factor1=asymp(1:nvar,i);
      sg1bd = sgpbida(1:nvar).*factor1;
      sg1bdinv = 1./sg1bd;
      %
      Hptd(1:nvar,1:nvar)=diag(sg1bd);
      Hptdinv(1:nvar,1:nvar)=diag(sg1bdinv);
      if lags>1
         factorpp=asymp(nvar+1:ncoef-1,i);
         sgpp_cbd = sgppbd(1:ncoef-nvar-1) .* factorpp;
         sgpp_cbdinv = 1./sgpp_cbd;
         Hptd(nvar+1:ncoef-1,nvar+1:ncoef-1)=diag(sgpp_cbd);
         Hptdinv(nvar+1:ncoef-1,nvar+1:ncoef-1)=diag(sgpp_cbdinv);
               % condtional on A0i, H_plus_tilde
      end
   end
   %---------------
   % The dummy observation prior affects only the prior covariance of A+|A0,
   %    but not the covariance of A0.  See pp.69a-69b for the proof.
   %---------------
   if indxDummy   % Dummy observations as part of the explicit prior.
      Hptdinv2 = Hptdinv + xtxbar;  % Rename Hptdinv to Hptdinv2 because we want to keep Hptdinv diagonal in the next loop of i.
   else
      Hptdinv2 = Hptdinv;
   end
   if (indxScaleStates)
      Hptldcell_inv{i}=(Vi{i}'*kron(eye(nStates),Hptdinv2/nStates))*Vi{i};
   else
      Hptldcell_inv{i}=(Vi{i}'*kron(eye(nStates),Hptdinv2))*Vi{i};
   end
   %Hptdinv_3 = kron(eye(nStates),Hptdinv);   % ?????
end


