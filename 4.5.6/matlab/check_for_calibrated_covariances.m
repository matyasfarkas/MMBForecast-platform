function estim_params=check_for_calibrated_covariances(xparam1,estim_params,M)
% function check_for_calibrated_covariances(xparam1,estim_params,M)
% find calibrated covariances to consider during estimation
% Inputs
%   -xparam1        [vector] parameters to be estimated
%   -estim_params   [structure] describing parameters to be estimated
%   -M              [structure] describing the model
%
% Outputs
%   -estim_params   [structure] describing parameters to be estimated
%
% Notes: M is local to this function and not updated when calling
% set_all_parameters

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
Sigma_e_calibrated=M.Sigma_e;
H_calibrated=M.H;
%check covariance for structural errors
covariance_pos=find(tril(Sigma_e_calibrated,-1)); %off-diagonal elements set by covariances before updating correlation matrix to reflect estimated covariances
covariance_pos_ME=find(tril(H_calibrated,-1)); %off-diagonal elements set by covariances before updating correlation matrix to reflect estimated covariances

%locally updated M
M = set_all_parameters(xparam1,estim_params,M);

correlation_pos=find(tril(M.Correlation_matrix,-1)); %off-diagonal elements set by correlations after accounting for estimation
calibrated_covariance_pos=covariance_pos(~ismember(covariance_pos,correlation_pos));
if any(calibrated_covariance_pos)
    [rows, columns]=ind2sub(size(M.Sigma_e),calibrated_covariance_pos); %find linear indices of lower triangular covariance entries
    estim_params.calibrated_covariances.position=[calibrated_covariance_pos;sub2ind(size(M.Sigma_e),columns,rows)]; %get linear entries of upper triangular parts
    estim_params.calibrated_covariances.cov_value=Sigma_e_calibrated(estim_params.calibrated_covariances.position);
end

correlation_pos_ME=find(tril(M.Correlation_matrix_ME,-1)); %off-diagonal elements set by correlations after accounting for estimation
calibrated_covariance_pos_ME=covariance_pos_ME(~ismember(covariance_pos_ME,correlation_pos_ME));
if any(calibrated_covariance_pos_ME)
    [rows, columns]=ind2sub(size(M.H),calibrated_covariance_pos_ME); %find linear indices of lower triangular covariance entries
    estim_params.calibrated_covariances_ME.position=[calibrated_covariance_pos_ME;sub2ind(size(M.H),columns,rows)]; %get linear entries of upper triangular parts
    estim_params.calibrated_covariances_ME.cov_value=H_calibrated(estim_params.calibrated_covariances_ME.position);
end
