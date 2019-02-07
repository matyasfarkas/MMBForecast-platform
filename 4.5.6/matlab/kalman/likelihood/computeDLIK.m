function [Da,DP,DLIK,D2a,D2P,Hesst] = computeDLIK(k,tmp,Z,Zflag,v,T,K,P,iF,Da,DYss,DT,DOm,DP,DH,notsteady,D2a,D2Yss,D2T,D2Om,D2P)

% Copyright (C) 2012-2017 Dynare Team
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

% AUTHOR(S) marco.ratto@jrc.ec.europa.eu

persistent DK DF D2K D2F

if notsteady
    if Zflag
        [DK,DF,DP1] = computeDKalmanZ(T,DT,DOm,P,DP,DH,Z,iF,K);
        if nargout>4
            [D2K,D2F,D2P] = computeD2KalmanZ(T,DT,D2T,D2Om,P,DP,D2P,DH,Z,iF,K,DK);
        end
    else
        [DK,DF,DP1] = computeDKalman(T,DT,DOm,P,DP,DH,Z,iF,K);
        if nargout>4
            [D2K,D2F,D2P] = computeD2Kalman(T,DT,D2T,D2Om,P,DP,D2P,DH,Z,iF,K,DK);
        end
    end
    DP=DP1;
    clear DP1
else
    DP=DP;
    if nargout>4
        D2P=D2P;
    end
end

Dv=zeros(length(v),k);
% D2v=zeros(length(v),k,k);
for ii = 1:k
    if Zflag
        Dv(:,ii)   = -Z*Da(:,ii) - Z*DYss(:,ii);
        %         if nargout>4,
        %             for jj = 1:ii
        %                 D2v(:,jj,ii)  = -Z*D2Yss(:,jj,ii)  - Z*D2a(:,jj,ii);
        %                 D2v(:,ii,jj) = D2v(:,jj,ii);
        %             end
        %         end
    else
        Dv(:,ii)   = -Da(Z,ii) - DYss(Z,ii);
        %         if nargout>4,
        %             for jj = 1:ii
        %                 D2v(:,jj,ii)  = -D2Yss(Z,jj,ii)  - D2a(Z,jj,ii);
        %                 D2v(:,ii,jj) = D2v(:,jj,ii);
        %             end
        %         end
    end
end

Hesst = zeros(k,k);
DLIK=zeros(k,1);
jcount=0;
for ii = 1:k
    %  dai = da(:,:,ii);
    dKi  = DK(:,:,ii);
    dtmp(:,ii) = Da(:,ii)+dKi*v+K*Dv(:,ii);

    if nargout>4
        diFi = -iF*DF(:,:,ii)*iF;
        for jj = 1:ii
            jcount=jcount+1;
            dFj    = DF(:,:,jj);
            diFj   = -iF*DF(:,:,jj)*iF;
            dKj  = DK(:,:,jj);
            d2Kij  = D2K(:,:,jj,ii);
            d2Fij  = D2F(:,:,jj,ii);
            d2iFij = -diFi*dFj*iF -iF*d2Fij*iF -iF*dFj*diFi;
            %             dtmpj = Da(:,jj)+dKj*v+K*Dv(:,jj);

            %             d2vij  = D2v(:,ii,jj);
            if Zflag
                d2vij  = -Z*D2Yss(:,jj,ii)  - Z*D2a(:,jj,ii);
            else
                d2vij  = -D2Yss(Z,jj,ii)  - D2a(Z,jj,ii);
            end
            d2tmpij = D2a(:,jj,ii) + d2Kij*v + dKj*Dv(:,ii) + dKi*Dv(:,jj) + K*d2vij;
            D2a(:,jj,ii) = reshape(D2T(:,jcount),size(T))*tmp + DT(:,:,jj)*dtmp(:,ii) + DT(:,:,ii)*dtmp(:,jj) + T*d2tmpij;
            D2a(:,ii,jj) = D2a(:,jj,ii);

            if nargout==6
                Hesst(ii,jj) = getHesst_ij(v,Dv(:,ii),Dv(:,jj),d2vij,iF,diFi,diFj,d2iFij,dFj,d2Fij);
            end
        end
    end

    Da(:,ii)   = DT(:,:,ii)*tmp + T*dtmp(:,ii);
    DLIK(ii,1)  = trace( iF*DF(:,:,ii) ) + 2*Dv(:,ii)'*iF*v - v'*(iF*DF(:,:,ii)*iF)*v;
end

if nargout==4
    %         Hesst(ii,jj) = getHesst_ij(v,Dv(:,ii),Dv(:,jj),0,iF,diFi,diFj,0,dFj,0);
    vecDPmf = reshape(DF,[],k);
    D2a = 2*Dv'*iF*Dv + (vecDPmf' * kron(iF,iF) * vecDPmf);
    %     for ii = 1:k
    %
    %         diFi = -iF*DF(:,:,ii)*iF;
    %         for jj = 1:ii
    %             dFj    = DF(:,:,jj);
    %             diFj   = -iF*DF(:,:,jj)*iF;
    %
    %             Hesst(ii,jj) = getHesst_ij(v*0,Dv(:,ii),Dv(:,jj),v*0,iF,diFi,diFj,0,-dFj,0);
    %         end
    %     end
end

% end of computeDLIK

function Hesst_ij = getHesst_ij(e,dei,dej,d2eij,iS,diSi,diSj,d2iSij,dSj,d2Sij)
% computes (i,j) term in the Hessian

Hesst_ij = trace(diSi*dSj + iS*d2Sij) + e'*d2iSij*e + 2*(dei'*diSj*e + dei'*iS*dej + e'*diSi*dej + e'*iS*d2eij);

% end of getHesst_ij

function [DK,DF,DP1] = computeDKalman(T,DT,DOm,P,DP1,DH,Z,iF,K)

k      = size(DT,3);
tmp    = P-K*P(Z,:);
DF = zeros([size(iF),k]);
DK = zeros([size(K),k]);
% DP1 = zeros([size(P),k]);

for ii = 1:k
    DF(:,:,ii)  = DP1(Z,Z,ii) + DH(:,:,ii);
    DiF = -iF*DF(:,:,ii)*iF;
    DK(:,:,ii)  = DP1(:,Z,ii)*iF + P(:,Z)*DiF;
    Dtmp        = DP1(:,:,ii) - DK(:,:,ii)*P(Z,:) - K*DP1(Z,:,ii);
    DP1(:,:,ii) = DT(:,:,ii)*tmp*T' + T*Dtmp*T' + T*tmp*DT(:,:,ii)' + DOm(:,:,ii);
end

% end of computeDKalman

function [DK,DF,DP1] = computeDKalmanZ(T,DT,DOm,P,DP,DH,Z,iF,K)

k      = size(DT,3);
tmp    = P-K*Z*P;
DF = zeros([size(iF),k]);
DK = zeros([size(K),k]);
DP1 = zeros([size(P),k]);

for ii = 1:k
    DF(:,:,ii)  = Z*DP(:,:,ii)*Z + DH(:,:,ii);
    DiF = -iF*DF(:,:,ii)*iF;
    DK(:,:,ii)  = DP(:,:,ii)*Z*iF + P(:,:)*Z*DiF;
    Dtmp        = DP(:,:,ii) - DK(:,:,ii)*Z*P(:,:) - K*Z*DP(:,:,ii);
    DP1(:,:,ii) = DT(:,:,ii)*tmp*T' + T*Dtmp*T' + T*tmp*DT(:,:,ii)' + DOm(:,:,ii);
end

% end of computeDKalmanZ

function [d2K,d2S,d2P1] = computeD2Kalman(A,dA,d2A,d2Om,P0,dP0,d2P1,DH,Z,iF,K0,dK0)
% computes the second derivatives of the Kalman matrices
% note: A=T in main func.

k      = size(dA,3);
tmp    = P0-K0*P0(Z,:);
[ns,no] = size(K0);

% CPC = C*P0*C'; CPC = .5*(CPC+CPC');iF = inv(CPC);
% APC = A*P0*C';
% APA = A*P0*A';


d2K  = zeros(ns,no,k,k);
d2S  = zeros(no,no,k,k);
% d2P1 = zeros(ns,ns,k,k);

jcount=0;
for ii = 1:k
    dAi = dA(:,:,ii);
    dFi = dP0(Z,Z,ii);
    %     d2Omi = d2Om(:,:,ii);
    diFi = -iF*dFi*iF;
    dKi = dK0(:,:,ii);
    for jj = 1:ii
        jcount=jcount+1;
        dAj = dA(:,:,jj);
        dFj = dP0(Z,Z,jj);
        %         d2Omj = d2Om(:,:,jj);
        dFj = dP0(Z,Z,jj);
        diFj = -iF*dFj*iF;
        dKj = dK0(:,:,jj);

        d2Aij = reshape(d2A(:,jcount),[ns ns]);
        d2Pij = dyn_unvech(d2P1(:,jcount));
        d2Omij = dyn_unvech(d2Om(:,jcount));

        % second order

        d2Fij = d2Pij(Z,Z) ;

        %     d2APC = d2Aij*P0*C' + A*d2Pij*C' + A*P0*d2Cij' + dAi*dPj*C' + dAj*dPi*C' + A*dPj*dCi' + A*dPi*dCj' + dAi*P0*dCj' + dAj*P0*dCi';
        d2APC = d2Pij(:,Z);

        d2iF = -diFi*dFj*iF -iF*d2Fij*iF -iF*dFj*diFi;

        d2Kij= d2Pij(:,Z)*iF + P0(:,Z)*d2iF + dP0(:,Z,jj)*diFi + dP0(:,Z,ii)*diFj;

        d2KCP = d2Kij*P0(Z,:) + K0*d2Pij(Z,:) + dKi*dP0(Z,:,jj) + dKj*dP0(Z,:,ii) ;

        dtmpi        = dP0(:,:,ii) - dK0(:,:,ii)*P0(Z,:) - K0*dP0(Z,:,ii);
        dtmpj        = dP0(:,:,jj) - dK0(:,:,jj)*P0(Z,:) - K0*dP0(Z,:,jj);
        d2tmp = d2Pij - d2KCP;

        d2AtmpA = d2Aij*tmp*A' + A*d2tmp*A' + A*tmp*d2Aij' + dAi*dtmpj*A' + dAj*dtmpi*A' + A*dtmpj*dAi' + A*dtmpi*dAj' + dAi*tmp*dAj' + dAj*tmp*dAi';

        d2K(:,:,ii,jj)  = d2Kij; %#ok<NASGU>
        d2P1(:,jcount) = dyn_vech(d2AtmpA  + d2Omij);  %#ok<*NASGU>
        d2S(:,:,ii,jj)  = d2Fij;
        d2K(:,:,jj,ii)  = d2Kij; %#ok<NASGU>
                                 %     d2P1(:,:,jj,ii) = d2AtmpA  + d2Omij;  %#ok<*NASGU>
        d2S(:,:,jj,ii)  = d2Fij;
        %     d2iS(:,:,ii,jj) = d2iF;
    end
end

% end of computeD2Kalman

function [d2K,d2S,d2P1] = computeD2KalmanZ(A,dA,d2A,d2Om,P0,dP0,d2P1,DH,Z,iF,K0,dK0)
% computes the second derivatives of the Kalman matrices
% note: A=T in main func.

k      = size(dA,3);
tmp    = P0-K0*Z*P0(:,:);
[ns,no] = size(K0);

% CPC = C*P0*C'; CPC = .5*(CPC+CPC');iF = inv(CPC);
% APC = A*P0*C';
% APA = A*P0*A';


d2K  = zeros(ns,no,k,k);
d2S  = zeros(no,no,k,k);
% d2P1 = zeros(ns,ns,k,k);

jcount=0;
for ii = 1:k
    dAi = dA(:,:,ii);
    dFi = Z*dP0(:,:,ii)*Z;
    %     d2Omi = d2Om(:,:,ii);
    diFi = -iF*dFi*iF;
    dKi = dK0(:,:,ii);
    for jj = 1:ii
        jcount=jcount+1;
        dAj = dA(:,:,jj);
        dFj = Z*dP0(:,:,jj)*Z;
        %         d2Omj = d2Om(:,:,jj);
        dFj = Z*dP0(:,:,jj)*Z;
        diFj = -iF*dFj*iF;
        dKj = dK0(:,:,jj);

        d2Aij = reshape(d2A(:,jcount),[ns ns]);
        d2Pij = dyn_unvech(d2P1(:,jcount));
        d2Omij = dyn_unvech(d2Om(:,jcount));

        % second order

        d2Fij = Z*d2Pij(:,:)*Z ;

        %     d2APC = d2Aij*P0*C' + A*d2Pij*C' + A*P0*d2Cij' + dAi*dPj*C' + dAj*dPi*C' + A*dPj*dCi' + A*dPi*dCj' + dAi*P0*dCj' + dAj*P0*dCi';
        d2APC = d2Pij(:,:)*Z;

        d2iF = -diFi*dFj*iF -iF*d2Fij*iF -iF*dFj*diFi;

        d2Kij= d2Pij(:,:)*Z*iF + P0(:,:)*Z*d2iF + dP0(:,:,jj)*Z*diFi + dP0(:,:,ii)*Z*diFj;

        d2KCP = d2Kij*Z*P0(:,:) + K0*Z*d2Pij(:,:) + dKi*Z*dP0(:,:,jj) + dKj*Z*dP0(:,:,ii) ;

        dtmpi        = dP0(:,:,ii) - dK0(:,:,ii)*Z*P0(:,:) - K0*Z*dP0(:,:,ii);
        dtmpj        = dP0(:,:,jj) - dK0(:,:,jj)*Z*P0(:,:) - K0*Z*dP0(:,:,jj);
        d2tmp = d2Pij - d2KCP;

        d2AtmpA = d2Aij*tmp*A' + A*d2tmp*A' + A*tmp*d2Aij' + dAi*dtmpj*A' + dAj*dtmpi*A' + A*dtmpj*dAi' + A*dtmpi*dAj' + dAi*tmp*dAj' + dAj*tmp*dAi';

        d2K(:,:,ii,jj)  = d2Kij; %#ok<NASGU>
        d2P1(:,jcount) = dyn_vech(d2AtmpA  + d2Omij);  %#ok<*NASGU>
        d2S(:,:,ii,jj)  = d2Fij;
        %     d2iS(:,:,ii,jj) = d2iF;
        d2K(:,:,jj,ii)  = d2Kij; %#ok<NASGU>
                                 %     d2P1(:,:,jj,ii) = d2AtmpA  + d2Omij;  %#ok<*NASGU>
        d2S(:,:,jj,ii)  = d2Fij;
    end
end

% end of computeD2KalmanZ
