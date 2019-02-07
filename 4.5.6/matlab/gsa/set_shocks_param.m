function set_shocks_param(xparam1)
% function set_shocks_param(xparam1)
% Set the structural and measurement error variances and covariances

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

global estim_params_ M_

nvx = estim_params_.nvx;
ncx = estim_params_.ncx;
nvn = estim_params_.nvn;
ncn = estim_params_.ncn;

Sigma_e = M_.Sigma_e;
Correlation_matrix = M_.Correlation_matrix;

H = M_.H;
Correlation_matrix_ME = M_.Correlation_matrix_ME;
% setting shocks variance on the diagonal of Covariance matrix; used later
% for updating covariances
if nvx
    var_exo = estim_params_.var_exo;
    for i=1:nvx
        k =var_exo(i,1);
        Sigma_e(k,k) = xparam1(i)^2;
    end
end
% update offset
offset = nvx;

% setting measument error variance; on the diagonal of Covariance matrix; used later
% for updating covariances
if nvn
    for i=1:nvn
        k = estim_params_.nvn_observable_correspondence(i,1);
        H(k,k) = xparam1(i+offset)^2;
    end
end

% update offset
offset = nvx+nvn;

% setting shocks covariances
if ncx
    corrx = estim_params_.corrx;
    for i=1:ncx
        k1 = corrx(i,1);
        k2 = corrx(i,2);
        Correlation_matrix(k1,k2) = xparam1(i+offset);
        Correlation_matrix(k2,k1) = Correlation_matrix(k1,k2);
    end
end
%build covariance matrix from correlation matrix and variances already on
%diagonal
Sigma_e = diag(sqrt(diag(Sigma_e)))*Correlation_matrix*diag(sqrt(diag(Sigma_e)));
%if calibrated covariances, set them now to their stored value
if isfield(estim_params_,'calibrated_covariances')
    Sigma_e(estim_params_.calibrated_covariances.position)=estim_params_.calibrated_covariances.cov_value;
end
% update offset
offset = nvx+nvn+ncx;

% setting measurement error covariances
if ncn
    corrn_observable_correspondence = estim_params_.corrn_observable_correspondence;
    for i=1:ncn
        k1 = corrn_observable_correspondence(i,1);
        k2 = corrn_observable_correspondence(i,2);
        Correlation_matrix_ME(k1,k2) = xparam1(i+offset);
        Correlation_matrix_ME(k2,k1) = Correlation_matrix_ME(k1,k2);
    end
end
%build covariance matrix from correlation matrix and variances already on
%diagonal
H = diag(sqrt(diag(H)))*Correlation_matrix_ME*diag(sqrt(diag(H)));
%if calibrated covariances, set them now to their stored value
if isfield(estim_params_,'calibrated_covariances_ME')
    H(estim_params_.calibrated_covariances_ME.position)=estim_params_.calibrated_covariances_ME.cov_value;
end


% updating matrices in M
if nvx || ncx
    M_.Sigma_e = Sigma_e;
    M_.Correlation_matrix=Correlation_matrix;
end
if nvn || ncn
    M_.H = H;
    M_.Correlation_matrix_ME=Correlation_matrix_ME;
end