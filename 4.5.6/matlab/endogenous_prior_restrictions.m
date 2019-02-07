function [info, info_irf, info_moment, data_irf, data_moment] = endogenous_prior_restrictions(T,R,Model,DynareOptions,DynareResults)
% Check for prior (sign) restrictions on irf's and theoretical moments
%
% INPUTS
%    T          [double]     n*n state space matrix
%    R          [double]     n*k matrix of shocks
%    Model      [structure]
%    DynareOptions [structure]
%    DynareResults [structure]

% OUTPUTS
%    info     [double]  check if prior restrictions are matched by the
%                       model and related info
%    info_irf [double] array of test checks for all individual irf restrictions
%    info_moment [double] array of test checks for all individual moment restrictions
%

% Copyright (C) 2013-2017 Dynare Team
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

info=[0 0];
info_irf=[];
info_moment=[];
data_irf=[];
data_moment=[];

endo_prior_restrictions.irf= DynareOptions.endogenous_prior_restrictions.irf;
endo_prior_restrictions.moment= DynareOptions.endogenous_prior_restrictions.moment;

if ~isempty(endo_prior_restrictions.irf)
    data_irf=cell(size(endo_prior_restrictions.irf,1),1);
    if DynareOptions.order>1
        error('The algorithm for prior (sign) restrictions on irf''s is currently restricted to first-order decision rules')
        return
    end
    varlist=Model.endo_names(DynareResults.dr.order_var,:);
    if isempty(T)
        [T,R,SteadyState,infox,Model,DynareOptions,DynareResults] = dynare_resolve(Model,DynareOptions,DynareResults);
    else % check if T and R are given in the restricted form!!!
        if size(T,1)<size(varlist,1)
            varlist=varlist(DynareResults.dr.restrict_var_list,:);
        end
        % check if endo_prior_restrictions.irf{:,1} variables are in varlist
        varlistok=1;
        for j=1:size(endo_prior_restrictions.irf,1)
            if isempty(strmatch(endo_prior_restrictions.irf{j,1},varlist,'exact'))
                varlistok=0;
            end
        end
        if ~varlistok
            varlist=Model.endo_names(DynareResults.dr.order_var,:);
            [T,R,SteadyState,infox,Model,DynareOptions,DynareResults] = dynare_resolve(Model,DynareOptions,DynareResults);
        end
    end
    NT=1;
    for j=1:size(endo_prior_restrictions.irf,1)
        NT=max(NT,max(endo_prior_restrictions.irf{j,3}));
    end
    info_irf=ones(size(endo_prior_restrictions.irf,1),2);
    for t=1:NT
        if ~DynareOptions.relative_irf
            RR = T^(t-1)*R*diag(sqrt(diag(Model.Sigma_e)));
        else
            RR = T^(t-1)*R*100;
        end
        for j=1:size(endo_prior_restrictions.irf,1)
            if endo_prior_restrictions.irf{j,3}~=t
                continue
            end
            iendo=strmatch(endo_prior_restrictions.irf{j,1},varlist,'exact');
            iexo=strmatch(endo_prior_restrictions.irf{j,2},Model.exo_names,'exact');
            data_irf{j}=[data_irf{j}; [t RR(iendo,iexo)]];
            if (RR(iendo,iexo)>endo_prior_restrictions.irf{j,4}(1)) && (RR(iendo,iexo)<endo_prior_restrictions.irf{j,4}(2))
                info_irf(j,:)=info_irf(j,:).*[0, 0];
            else
                if RR(iendo,iexo)<endo_prior_restrictions.irf{j,4}(1)
                    delt = (RR(iendo,iexo)-endo_prior_restrictions.irf{j,4}(1))^2;
                else
                    delt = (RR(iendo,iexo)-endo_prior_restrictions.irf{j,4}(2))^2;
                end
                info_irf(j,:)=info_irf(j,:).*[49, delt];
            end
        end
    end
    if any(info_irf)
        info=[49,sum(info_irf(:,2))];
    end

end

if ~isempty(endo_prior_restrictions.moment)
    if DynareOptions.order>1
        error('The algorithm for prior (sign) restrictions on moments is currently restricted to first-order decision rules')
        return
    end
    data_moment=cell(size(endo_prior_restrictions.moment,1),1);
    var_list_=endo_prior_restrictions.moment{1,1};
    for  j=1:size(endo_prior_restrictions.moment,1)
        tmp=endo_prior_restrictions.moment{j,1};
        if ~ismember(tmp,cellstr(var_list_))
            var_list_ = char(var_list_, tmp);
        end
        tmp=endo_prior_restrictions.moment{j,2};
        if ~ismember(tmp,cellstr(var_list_))
            var_list_ = char(var_list_, tmp);
        end
    end
    NTmax=0;
    NTmin=0;
    for j=1:size(endo_prior_restrictions.moment,1)
        NTmax=max(NTmax,max(endo_prior_restrictions.moment{j,3}));
        NTmin=min(NTmin,min(endo_prior_restrictions.moment{j,3}));
    end
    info_moment=ones(size(endo_prior_restrictions.moment,1),2);
    nvar = size(var_list_,1);
    ivar=zeros(nvar,1);
    for i=1:nvar
        i_tmp = strmatch(var_list_(i,:),Model.endo_names,'exact');
        if isempty(i_tmp)
            error (['One of the variable specified does not exist'])
        else
            ivar(i) = i_tmp;
        end
    end
    DynareOptions.ar = max(abs(NTmin),NTmax);
    [gamma_y,stationary_vars] = th_autocovariances(DynareResults.dr, ivar, Model, DynareOptions,1);
    for t=NTmin:NTmax
        RR = gamma_y{abs(t)+1};
        if t==0
            RR = RR./(sqrt(diag(RR))*sqrt(diag(RR))')-eye(nvar)+diag(diag(gamma_y{t+1})); % becomes correlation
        end
        for j=1:size(endo_prior_restrictions.moment,1)
            if endo_prior_restrictions.moment{j,3}~=t
                continue
            end
            iendo1 = strmatch(endo_prior_restrictions.moment{j,1},var_list_,'exact');
            iendo2 = strmatch(endo_prior_restrictions.moment{j,2},var_list_,'exact');
            if t>0
                tmp0 = iendo1;
                iendo1=iendo2;
                iendo2=tmp0;
            end
            data_moment{j}=[data_moment{j}; [t RR(iendo1,iendo2)]];
            if (RR(iendo1,iendo2)>endo_prior_restrictions.moment{j,4}(1)) && (RR(iendo1,iendo2)<endo_prior_restrictions.moment{j,4}(2))
                info_moment(j,:)=info_moment(j,:).*[0, 0];
            else
                if RR(iendo1,iendo2)<endo_prior_restrictions.moment{j,4}(1)
                    delt = (RR(iendo1,iendo2)-endo_prior_restrictions.moment{j,4}(1))^2;
                else
                    delt = (RR(iendo1,iendo2)-endo_prior_restrictions.moment{j,4}(2))^2;
                end
                info_moment(j,:)=info_moment(j,:).*[49, delt];
            end
        end
    end
    if any(info_moment)
        info=[49, info(2) + sum(info_moment(:,2))];
    end
end