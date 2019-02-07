function [condJ, ind0, indnoJ, ixnoJ, McoJ, PcoJ, jweak, jweak_pair] = identification_checks(JJ, hess_flag)
% function [condJ, ind0, indnoJ, ixnoJ, McoJ, PcoJ, jweak, jweak_pair] = identification_checks(JJ, hess_flag)
% checks for identification
%
% INPUTS
%    o JJ               [matrix] [output x nparams] IF hess_flag==0
%                                 derivatives of output w.r.t. parameters and shocks
%    o JJ               [matrix] [nparams x nparams] IF hess_flag==1
%                                 information matrix
%
% OUTPUTS
%    o cond             condition number of JJ
%    o ind0             [array] binary indicator for non-zero columns of H
%    o indnoJ           [matrix] index of non-identified params
%    o ixnoJ            number of rows in indnoJ
%    o Mco              [array] multicollinearity coefficients
%    o Pco              [matrix] pairwise correlations
%    o jweak            [binary array] gives 1 if the  parameter has Mco=1(with tolerance 1.e-10)
%    o jweak_pair       [binary matrix] gives 1 if a couple parameters has Pco=1(with tolerance 1.e-10)
%
% SPECIAL REQUIREMENTS
%    None

% Copyright (C) 2008-2017 Dynare Team
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

% My suggestion is to have the following steps for identification check in
% dynare:

% 1. check rank of JJ at theta
npar = size(JJ,2);
indnoJ = zeros(1,npar);

if size(JJ,1)>1
    ind1 = find(vnorm(JJ)>=eps); % take non-zero columns
else
    ind1 = find(abs(JJ)>=eps); % take non-zero columns
end
JJ1 = JJ(:,ind1);
[eu,ee2,ee1] = svd( JJ1, 0 );
condJ= cond(JJ1);
rankJ = rank(JJ);
rankJJ = rankJ;
icheck=0;
if npar>0 && (rankJ<npar)
    % search for singular values associated to ONE individual parameter
    ee0 = [rankJJ+1:length(ind1)];
    ind11=ones(length(ind1),1);
    for j=1:length(ee0)
        if length(find(abs(ee1(:,ee0(j))) > 1.e-3))==1
            icheck=1;
            ind11 = ind11.*(abs(ee1(:,ee0(j))) <= 1.e-3); % take non-zero columns
        end
    end
    ind1 = ind1(find(ind11)); % take non-zero columns
end

if icheck
    JJ1 = JJ(:,ind1);
    [eu,ee2,ee1] = svd( JJ1, 0 );
    condJ= cond(JJ1);
    rankJ = rank(JJ);
    rankJJ = rankJ;
end


% if hess_flag==0,
%     rankJJ = rank(JJ'*JJ);
% end

ind0 = zeros(1,npar);
ind0(ind1) = 1;

if hess_flag==0
    % find near linear dependence problems:
    McoJ = NaN(npar,1);
    for ii = 1:size(JJ1,2)
        McoJ(ind1(ii),:) = cosn([JJ1(:,ii),JJ1(:,find([1:1:size(JJ1,2)]~=ii))]);
    end
else
    deltaJ = sqrt(diag(JJ(ind1,ind1)));
    tildaJ = JJ(ind1,ind1)./((deltaJ)*(deltaJ'));
    McoJ(ind1,1)=(1-1./diag(inv(tildaJ)));
    rhoM=sqrt(1-McoJ);
    %     PcoJ=inv(tildaJ);
    PcoJ=NaN(npar,npar);
    PcoJ(ind1,ind1)=inv(JJ(ind1,ind1));
    sd=sqrt(diag(PcoJ));
    PcoJ = abs(PcoJ./((sd)*(sd')));
end


ixnoJ = 0;
if npar>0 && (rankJ<npar || rankJJ<npar || min(1-McoJ)<1.e-10)
    %         - find out which parameters are involved
    %   disp('Some parameters are NOT identified by the moments included in J')
    %   disp(' ')
    if length(ind1)<npar
        % parameters with zero column in JJ
        ixnoJ = ixnoJ + 1;
        indnoJ(ixnoJ,:) = (~ismember([1:npar],ind1));
    end
    ee0 = [rankJJ+1:length(ind1)];
    for j=1:length(ee0)
        % linearely dependent parameters in JJ
        ixnoJ = ixnoJ + 1;
        indnoJ(ixnoJ,ind1) = (abs(ee1(:,ee0(j))) > 1.e-3)';
    end
end

% here there is no exact linear dependence, but there are several
%     near-dependencies, mostly due to strong pairwise colliniearities, which can
%     be checked using

jweak=zeros(1,npar);
jweak_pair=zeros(npar,npar);

if hess_flag==0,
    PcoJ = NaN(npar,npar);

    for ii = 1:size(JJ1,2)
        PcoJ(ind1(ii),ind1(ii)) = 1;
        for jj = ii+1:size(JJ1,2)
            PcoJ(ind1(ii),ind1(jj)) = cosn([JJ1(:,ii),JJ1(:,jj)]);
            PcoJ(ind1(jj),ind1(ii)) = PcoJ(ind1(ii),ind1(jj));
        end
    end

    for j=1:npar
        if McoJ(j)>(1-1.e-10)
            jweak(j)=1;
            [ipair, jpair] = find(PcoJ(j,j+1:end)>(1-1.e-10));
            for jx=1:length(jpair)
                jweak_pair(j, jpair(jx)+j)=1;
                jweak_pair(jpair(jx)+j, j)=1;
            end
        end
    end
end

jweak_pair=dyn_vech(jweak_pair)';
