function oo_ = covariance_mc_analysis(NumberOfSimulations,type,dname,fname,vartan,nvar,var1,var2,mh_conf_sig,oo_,options_)
% This function analyses the (posterior or prior) distribution of the
% endogenous variables' covariance matrix.
%
% INPUTS
%   NumberOfSimulations     [integer]           scalar, number of simulations.
%   type                    [string]            'prior' or 'posterior'
%   dname                   [string]            directory name where to save
%   fname                   [string]            name of the mod-file
%   vartan                  [char]              array of characters (with nvar rows).
%   nvar                    [integer]           nvar is the number of stationary variables.
%   var1                    [string]            name of the first variable
%   var2                    [string]            name of the second variable
%   mh_conf_sig             [double]            2 by 1 vector with upper
%                                               and lower bound of HPD intervals
%   oo_                     [structure]         Dynare structure where the results are saved.
%   options_                [structure]         Dynare options structure
%
% OUTPUTS
%   oo_                     [structure]        Dynare structure where the results are saved.

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

if strcmpi(type,'posterior')
    TYPE = 'Posterior';
    PATH = [dname '/metropolis/'];
else
    TYPE = 'Prior';
    PATH = [dname '/prior/moments/'];
end

indx1 = check_name(vartan,var1);
if isempty(indx1)
    disp([ type '_analysis:: ' var1 ' is not a stationary endogenous variable!'])
    return
end
if ~isempty(var2)
    indx2 = check_name(vartan,var2);
    if isempty(indx2)
        disp([ type '_analysis:: ' var2 ' is not a stationary endogenous variable!'])
        return
    end
else
    indx2 = indx1;
    var2 = var1;
end

var1=deblank(var1);
var2=deblank(var2);

if isfield(oo_,[ TYPE 'TheoreticalMoments'])
    temporary_structure = oo_.([TYPE, 'TheoreticalMoments']);
    if isfield(temporary_structure,'dsge')
        temporary_structure = oo_.([TYPE, 'TheoreticalMoments']).dsge;
        if isfield(temporary_structure,'covariance')
            temporary_structure = oo_.([TYPE, 'TheoreticalMoments']).dsge.covariance.Mean;
            if isfield(temporary_structure,var1)
                temporary_structure_1 = oo_.([TYPE, 'TheoreticalMoments']).dsge.covariance.Mean.(var1);
                if isfield(temporary_structure_1,var2)
                    % Nothing to do (the covariance matrix is symmetric!).
                    return
                end
            else
                if isfield(temporary_structure,var2)
                    temporary_structure_2 = oo_.([TYPE, 'TheoreticalMoments']).dsge.covariance.Mean.(var2);
                    if isfield(temporary_structure_2,var1)
                        % Nothing to do (the covariance matrix is symmetric!).
                        return
                    end
                end
            end
        end
    end
end

ListOfFiles = dir([ PATH  fname '_' TYPE '2ndOrderMoments*.mat']);
i1 = 1; tmp = zeros(NumberOfSimulations,1);
if options_.contemporaneous_correlation
    tmp_corr_mat = zeros(NumberOfSimulations,1);
    cov_pos=symmetric_matrix_index(indx1,indx2,nvar);
    var_pos_1=symmetric_matrix_index(indx1,indx1,nvar);
    var_pos_2=symmetric_matrix_index(indx2,indx2,nvar);
end
for file = 1:length(ListOfFiles)
    load([ PATH ListOfFiles(file).name ]);
    i2 = i1 + rows(Covariance_matrix) - 1;
    tmp(i1:i2) = Covariance_matrix(:,symmetric_matrix_index(indx1,indx2,nvar));
    if options_.contemporaneous_correlation
        temp=Covariance_matrix(:,cov_pos)./(sqrt(Covariance_matrix(:,var_pos_1)).*sqrt(Covariance_matrix(:,var_pos_2)));
        temp(Covariance_matrix(:,cov_pos)==0)=0; %filter out 0 correlations that would result in 0/0
        tmp_corr_mat(i1:i2)=temp;
    end
    i1 = i2+1;
end

if options_.estimation.moments_posterior_density.indicator
    [p_mean, p_median, p_var, hpd_interval, p_deciles, density] = ...
        posterior_moments(tmp,1,mh_conf_sig);
    oo_.([TYPE, 'TheoreticalMoments']).dsge.covariance.density.(var1).(var2) = density;
else
    [p_mean, p_median, p_var, hpd_interval, p_deciles] = ...
        posterior_moments(tmp,0,mh_conf_sig);
end
oo_.([TYPE, 'TheoreticalMoments']).dsge.covariance.Mean.(var1).(var2) = p_mean;
oo_.([TYPE, 'TheoreticalMoments']).dsge.covariance.Median.(var1).(var2) = p_median;
oo_.([TYPE, 'TheoreticalMoments']).dsge.covariance.Variance.(var1).(var2) = p_var;
oo_.([TYPE, 'TheoreticalMoments']).dsge.covariance.HPDinf.(var1).(var2) = hpd_interval(1);
oo_.([TYPE, 'TheoreticalMoments']).dsge.covariance.HPDsup.(var1).(var2) = hpd_interval(2);
oo_.([TYPE, 'TheoreticalMoments']).dsge.covariance.deciles.(var1).(var2) = p_deciles;

if options_.contemporaneous_correlation
    if options_.estimation.moments_posterior_density.indicator
        [p_mean, p_median, p_var, hpd_interval, p_deciles, density] = ...
            posterior_moments(tmp_corr_mat,1,mh_conf_sig);
        oo_.([TYPE, 'TheoreticalMoments']).dsge.contemporeaneous_correlation.density.(var1).(var2) = density;
    else
        [p_mean, p_median, p_var, hpd_interval, p_deciles] = ...
            posterior_moments(tmp_corr_mat,0,mh_conf_sig);
    end
    oo_.([TYPE, 'TheoreticalMoments']).dsge.contemporeaneous_correlation.Mean.(var1).(var2) = p_mean;
    oo_.([TYPE, 'TheoreticalMoments']).dsge.contemporeaneous_correlation.Median.(var1).(var2) = p_median;
    oo_.([TYPE, 'TheoreticalMoments']).dsge.contemporeaneous_correlation.Variance.(var1).(var2) = p_var;
    oo_.([TYPE, 'TheoreticalMoments']).dsge.contemporeaneous_correlation.HPDinf.(var1).(var2) = hpd_interval(1);
    oo_.([TYPE, 'TheoreticalMoments']).dsge.contemporeaneous_correlation.HPDsup.(var1).(var2) = hpd_interval(2);
    oo_.([TYPE, 'TheoreticalMoments']).dsge.contemporeaneous_correlation.deciles.(var1).(var2) = p_deciles;
end
