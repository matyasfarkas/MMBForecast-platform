function [Hess] = get_Hessian(T,R,Q,H,P,Y,DT,DYss,DOm,DH,DP,D2T,D2Yss,D2Om,D2H,D2P,start,mf,kalman_tol,riccati_tol)
% function [Hess] = get_Hessian(T,R,Q,H,P,Y,DT,DYss,DOm,DH,DP,D2T,D2Yss,D2Om,D2H,D2P,start,mf,kalman_tol,riccati_tol)
%
% computes the hessian matrix of the log-likelihood function of
% a state space model (notation as in kalman_filter.m in DYNARE
% Thanks to  Nikolai Iskrev
%
% NOTE: the derivative matrices (DT,DR ...) are 3-dim. arrays with last
% dimension equal to the number of structural parameters
% NOTE: the derivative matrices (D2T,D2Om ...) are 4-dim. arrays with last
% two dimensions equal to the number of structural parameters

% Copyright (C) 2011-2017 Dynare Team
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


k = size(DT,3);                                 % number of structural parameters
smpl = size(Y,2);                               % Sample size.
pp   = size(Y,1);                               % Maximum number of observed variables.
mm   = size(T,2);                               % Number of state variables.
a    = zeros(mm,1);                             % State vector.
Om   = R*Q*transpose(R);                        % Variance of R times the vector of structural innovations.
t    = 0;                                       % Initialization of the time index.
oldK = 0;
notsteady   = 1;                                % Steady state flag.
F_singular  = 1;

Hess  = zeros(k,k);                             % Initialization of the Hessian
Da    = zeros(mm,k);                             % State vector.
Dv = zeros(length(mf),k);
D2a    = zeros(mm,k,k);                             % State vector.
D2v = zeros(length(mf),k,k);

C = zeros(length(mf),mm);
for ii=1:length(mf); C(ii,mf(ii))=1;end         % SELECTION MATRIX IN MEASUREMENT EQ. (FOR WHEN IT IS NOT CONSTANT)
dC = zeros(length(mf),mm,k);
d2C = zeros(length(mf),mm,k,k);

s   = zeros(pp,1);                      % CONSTANT TERM IN MEASUREMENT EQ. (FOR WHEN IT IS NOT CONSTANT)
ds  = zeros(pp,1,k);
d2s = zeros(pp,1,k,k);

%     for ii = 1:k
%         DOm = DR(:,:,ii)*Q*transpose(R) + R*DQ(:,:,ii)*transpose(R) + R*Q*transpose(DR(:,:,ii));
%     end

while notsteady & t<smpl
    t  = t+1;
    v  = Y(:,t)-a(mf);
    F  = P(mf,mf) + H;
    if rcond(F) < kalman_tol
        if ~all(abs(F(:))<kalman_tol)
            return
        else
            a = T*a;
            P = T*P*transpose(T)+Om;
        end
    else
        F_singular = 0;
        iF     = inv(F);
        K      = P(:,mf)*iF;

        [DK,DF,DP1] = computeDKalman(T,DT,DOm,P,DP,DH,mf,iF,K);
        [D2K,D2F,D2P1] = computeD2Kalman(T,DT,D2T,D2Om,P,DP,D2P,DH,mf,iF,K,DK);
        tmp = (a+K*v);

        for ii = 1:k
            Dv(:,ii)   = -Da(mf,ii) - DYss(mf,ii);
            %  dai = da(:,:,ii);
            dKi  = DK(:,:,ii);
            diFi = -iF*DF(:,:,ii)*iF;
            dtmpi = Da(:,ii)+dKi*v+K*Dv(:,ii);


            for jj = 1:ii
                dFj    = DF(:,:,jj);
                diFj   = -iF*DF(:,:,jj)*iF;
                dKj  = DK(:,:,jj);
                d2Kij  = D2K(:,:,jj,ii);
                d2Fij  = D2F(:,:,jj,ii);
                d2iFij = -diFi*dFj*iF -iF*d2Fij*iF -iF*dFj*diFi;
                dtmpj = Da(:,jj)+dKj*v+K*Dv(:,jj);

                d2vij  = -D2Yss(mf,jj,ii)  - D2a(mf,jj,ii);
                d2tmpij = D2a(:,jj,ii) + d2Kij*v + dKj*Dv(:,ii) + dKi*Dv(:,jj) + K*d2vij;
                D2a(:,jj,ii) = D2T(:,:,jj,ii)*tmp + DT(:,:,jj)*dtmpi + DT(:,:,ii)*dtmpj + T*d2tmpij;

                Hesst(ii,jj) = getHesst_ij(v,Dv(:,ii),Dv(:,jj),d2vij,iF,diFi,diFj,d2iFij,dFj,d2Fij);
            end
            Da(:,ii)   = DT(:,:,ii)*tmp + T*dtmpi;
        end
        %                     vecDPmf = reshape(DP(mf,mf,:),[],k);
        %                     iPmf = inv(P(mf,mf));
        if t>=start
            Hess = Hess + Hesst;
        end
        a      = T*(a+K*v);
        P      = T*(P-K*P(mf,:))*transpose(T)+Om;
        DP     = DP1;
        D2P     = D2P1;
    end
    notsteady = max(max(abs(K-oldK))) > riccati_tol;
    oldK = K;
end

if F_singular
    error('The variance of the forecast error remains singular until the end of the sample')
end


if t < smpl
    t0 = t+1;
    while t < smpl
        t = t+1;
        v = Y(:,t)-a(mf);
        tmp = (a+K*v);
        for ii = 1:k
            Dv(:,ii)   = -Da(mf,ii)-DYss(mf,ii);
            dKi  = DK(:,:,ii);
            diFi = -iF*DF(:,:,ii)*iF;
            dtmpi = Da(:,ii)+dKi*v+K*Dv(:,ii);

            for jj = 1:ii
                dFj    = DF(:,:,jj);
                diFj   = -iF*DF(:,:,jj)*iF;
                dKj  = DK(:,:,jj);
                d2Kij  = D2K(:,:,jj,ii);
                d2Fij  = D2F(:,:,jj,ii);
                d2iFij = -diFi*dFj*iF -iF*d2Fij*iF -iF*dFj*diFi;
                dtmpj = Da(:,jj)+dKj*v+K*Dv(:,jj);

                d2vij  = -D2Yss(mf,jj,ii)  - D2a(mf,jj,ii);
                d2tmpij = D2a(:,jj,ii) + d2Kij*v + dKj*Dv(:,ii) + dKi*Dv(:,jj) + K*d2vij;
                D2a(:,jj,ii) = D2T(:,:,jj,ii)*tmp + DT(:,:,jj)*dtmpi + DT(:,:,ii)*dtmpj + T*d2tmpij;

                Hesst(ii,jj) = getHesst_ij(v,Dv(:,ii),Dv(:,jj),d2vij,iF,diFi,diFj,d2iFij,dFj,d2Fij);
            end
            Da(:,ii)   = DT(:,:,ii)*tmp + T*dtmpi;
        end
        if t>=start
            Hess = Hess + Hesst;
        end
        a = T*(a+K*v);
    end
    %         Hess = Hess + .5*(smpl+t0-1)*(vecDPmf' * kron(iPmf,iPmf) * vecDPmf);
    %         for ii = 1:k;
    %             for jj = 1:ii
    %              H(ii,jj) = trace(iPmf*(.5*DP(mf,mf,ii)*iPmf*DP(mf,mf,jj) + Dv(:,ii)*Dv(:,jj)'));
    %             end
    %         end
end

Hess = Hess + tril(Hess,-1)';

Hess = -Hess/2;
% end of main function

function Hesst_ij = getHesst_ij(e,dei,dej,d2eij,iS,diSi,diSj,d2iSij,dSj,d2Sij);
% computes (i,j) term in the Hessian

Hesst_ij = trace(diSi*dSj + iS*d2Sij) + e'*d2iSij*e + 2*(dei'*diSj*e + dei'*iS*dej + e'*diSi*dej + e'*iS*d2eij);

% end of getHesst_ij

function [DK,DF,DP1] = computeDKalman(T,DT,DOm,P,DP,DH,mf,iF,K)

k      = size(DT,3);
tmp    = P-K*P(mf,:);

for ii = 1:k
    DF(:,:,ii)  = DP(mf,mf,ii) + DH(:,:,ii);
    DiF(:,:,ii) = -iF*DF(:,:,ii)*iF;
    DK(:,:,ii)  = DP(:,mf,ii)*iF + P(:,mf)*DiF(:,:,ii);
    Dtmp        = DP(:,:,ii) - DK(:,:,ii)*P(mf,:) - K*DP(mf,:,ii);
    DP1(:,:,ii) = DT(:,:,ii)*tmp*T' + T*Dtmp*T' + T*tmp*DT(:,:,ii)' + DOm(:,:,ii);
end

% end of computeDKalman

function [d2K,d2S,d2P1] = computeD2Kalman(A,dA,d2A,d2Om,P0,dP0,d2P0,DH,mf,iF,K0,dK0)
% computes the second derivatives of the Kalman matrices
% note: A=T in main func.

k      = size(dA,3);
tmp    = P0-K0*P0(mf,:);
[ns,no] = size(K0);

% CPC = C*P0*C'; CPC = .5*(CPC+CPC');iF = inv(CPC);
% APC = A*P0*C';
% APA = A*P0*A';


d2K  = zeros(ns,no,k,k);
d2S  = zeros(no,no,k,k);
d2P1 = zeros(ns,ns,k,k);

for ii = 1:k
    dAi = dA(:,:,ii);
    dFi = dP0(mf,mf,ii);
    d2Omi = d2Om(:,:,ii);
    diFi = -iF*dFi*iF;
    dKi = dK0(:,:,ii);
    for jj = 1:k
        dAj = dA(:,:,jj);
        dFj = dP0(mf,mf,jj);
        d2Omj = d2Om(:,:,jj);
        dFj = dP0(mf,mf,jj);
        diFj = -iF*dFj*iF;
        dKj = dK0(:,:,jj);

        d2Aij = d2A(:,:,jj,ii);
        d2Pij = d2P0(:,:,jj,ii);
        d2Omij = d2Om(:,:,jj,ii);

        % second order

        d2Fij = d2Pij(mf,mf) ;

        %     d2APC = d2Aij*P0*C' + A*d2Pij*C' + A*P0*d2Cij' + dAi*dPj*C' + dAj*dPi*C' + A*dPj*dCi' + A*dPi*dCj' + dAi*P0*dCj' + dAj*P0*dCi';
        d2APC = d2Pij(:,mf);

        d2iF = -diFi*dFj*iF -iF*d2Fij*iF -iF*dFj*diFi;

        d2Kij= d2Pij(:,mf)*iF + P0(:,mf)*d2iF + dP0(:,mf,jj)*diFi + dP0(:,mf,ii)*diFj;

        d2KCP = d2Kij*P0(mf,:) + K0*d2Pij(mf,:) + dKi*dP0(mf,:,jj) + dKj*dP0(mf,:,ii) ;

        dtmpi        = dP0(:,:,ii) - dK0(:,:,ii)*P0(mf,:) - K0*dP0(mf,:,ii);
        dtmpj        = dP0(:,:,jj) - dK0(:,:,jj)*P0(mf,:) - K0*dP0(mf,:,jj);
        d2tmp = d2Pij - d2KCP;

        d2AtmpA = d2Aij*tmp*A' + A*d2tmp*A' + A*tmp*d2Aij' + dAi*dtmpj*A' + dAj*dtmpi*A' + A*dtmpj*dAi' + A*dtmpi*dAj' + dAi*tmp*dAj' + dAj*tmp*dAi';

        d2K(:,:,ii,jj)  = d2Kij; %#ok<NASGU>
        d2P1(:,:,ii,jj) = d2AtmpA  + d2Omij;  %#ok<*NASGU>
        d2S(:,:,ii,jj)  = d2Fij;
        %     d2iS(:,:,ii,jj) = d2iF;
    end
end

% end of computeD2Kalman
