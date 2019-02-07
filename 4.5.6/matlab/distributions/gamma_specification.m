function [a, b] = gamma_specification(mu, sigma2, lb, name)   % --*-- Unitary tests --*--

% Returns the hyperparameters of the gamma distribution given the expectation and variance.
%
% INPUTS
% - mu     [double]   Expectation of the Gamma random variable.
% - sigma2 [double]   Variance of the Gamma random variable.
% - lb     [double]   Lower bound of the domain (default is zero).
% - name   [string]   Name of the parameter (or random variable).
%
% OUTPUTS
% - a      [double]   First hyperparameter of the Gamma density (shape).
% - b      [double]   Second hyperparameter of the Gamma density (scale).

% Copyright (C) 2015-2017 Dynare Team
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

if nargin<3
    lb = 0;
end

if nargin>3 && ~isempty(name)
    name1 = sprintf('for %s ', name);
    name2 = sprintf(' (for %s)', name);
else
    name1 = '';
    name2 = '';
end

if mu<lb
    error('The prior expectation (%f) %scannot be smaller than the lower bound of the Gamma distribution (%f)!', mu, name1, lb)
end

if isinf(sigma2)
    error('The variance of the Gamma distribution has to be finite%s!', name2)
end

mu = mu-lb;
b  = sigma2/mu;
a  = mu/b;

%@test:1
%$ try
%$    [a, b] = gamma_specification(.5, 1);
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    t(2) = abs(0.5-a*b)<1e-12;
%$    t(3) = abs(1.0-a*b*b)<1e-12;
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ try
%$    [a, b] = gamma_specification(2, 1);
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    t(2) = abs(2.0-a*b)<1e-12;
%$    t(3) = abs(1.0-a*b*b)<1e-12;
%$ end
%$ T = all(t);
%@eof:2

%@test:3
%$ try
%$    [a, b] = gamma_specification(2, 1, 3);
%$    t(1) = false;
%$ catch
%$    t(1) = true;
%$ end
%$
%$ T = all(t);
%@eof:3

%@test:4
%$ try
%$    [a, b] = gamma_specification(2, Inf, 3);
%$    t(1) = false;
%$ catch
%$    t(1) = true;
%$ end
%$
%$ T = all(t);
%@eof:4