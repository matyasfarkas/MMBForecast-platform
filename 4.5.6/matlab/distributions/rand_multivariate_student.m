function draw = rand_multivariate_student(Mean,Sigma_upper_chol,df)
% function draw = rand_multivariate_student(Mean,Sigma_upper_chol,df)
% Pseudo random draws from a multivariate student distribution,
% with expectation Mean, variance Sigma*df/(df-2) and degrees of freedom df>0.
%
% INPUTS
%
%    Mean               [double]    1*n vector, expectation of the multivariate random variable.
%    Sigma_upper_chol   [double]    n*n matrix, upper triangular Cholesky decomposition of Sigma
%                                   (the covariance matrix up to a factor df/(df-2)).
%    df                 [integer]   degrees of freedom.
%
% OUTPUTS
%    draw               [double]    1*n vector drawn from a multivariate normal distribution with expectation Mean and
%                                   covariance Sigma.
%
%
% NOTE See Zellner (appendix B.2, 1971) for a definition.
%    Computes the t-distributed random numbers from
%       X = \mu + Y\sqrt{\frac{\nu}{U}}
%   where
%       Y~N(0,Sigma) with Sigma=Sigma_upper_chol'*Sigma_upper_chol
%       U~\Chi^2_{\nu}
%   The latter is constructed as the sum of \nu standard normals.

% Copyright (C) 2003-2017 Dynare Team
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

n = length(Mean);
draw = Mean + randn(1,n) * Sigma_upper_chol * sqrt(df/sum(randn(df,1).^2));
