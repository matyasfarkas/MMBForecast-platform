function [A,B] = kalman_transition_matrix(dr,iv,ic,exo_nbr)
%function [A,B] = kalman_transition_matrix(dr,iv,ic,exo_nbr)
% Builds the transition equation of the state space representation out of ghx and ghu for Kalman filter
%
% INPUTS
%    dr:      structure of decisions rules for stochastic simulations
%    iv:      selected variables
%    ic:      state variables position in the transition matrix columns
%    exo_nbr: number of exogenous variables
%
% OUTPUTS
%    A:       matrix of predetermined variables effects in linear solution (ghx)
%    B:       matrix of shocks effects in linear solution (ghu)
%
% SPECIAL REQUIREMENTS
%    none

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

n_iv = length(iv);

A = zeros(n_iv,n_iv);

A(:,ic) = dr.ghx(iv,:);

if nargout>1
    B = dr.ghu(iv,:);
end
