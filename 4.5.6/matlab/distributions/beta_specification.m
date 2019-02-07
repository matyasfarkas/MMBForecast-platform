function [a, b] = beta_specification(mu, sigma2, lb, ub, name)   % --*-- Unitary tests --*--

% Returns the hyperparameters of the beta distribution given the expectation and variance.
%
% INPUTS
% - mu     [double]   Expectation of the Gamma random variable.
% - sigma2 [double]   Variance of the Gamma random variable.
% - lb     [double]   Lower bound of the domain (default is zero).
% - ub     [double]   Upper bound of the domain (default is one).
%
% OUTPUTS
% - a      [double]   First hyperparameter of the Beta density.
% - b      [double]   Second hyperparameter of the Beta density.

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
    ub = 1;
end

if nargin<4
    ub = 1;
end

if nargin>4 && ~isempty(name)
    name1  = sprintf('of %s ', name);
    name2 =  sprintf(' (for %s)', name);
else
    name1 = '';
    name2 = '';
end

if mu<lb
    error(['The prior expectation (%f) %scannot be smaller than the lower bound of the Beta distribution (%f)!'], mu, name1, lb)
end

if mu>ub
    error('The prior expectation (%f) %scannot be greater than the upper bound of the Beta distribution (%f)!', mu, name1, ub)
end

len = ub-lb;

mu = (mu-lb)/len;
sigma2 = sigma2/(len*len);

if sigma2>(1-mu)*mu
    error('Beta prior%s. Given the declared prior expectation, prior lower and upper bounds, the prior std. has to be smaller than %f.',name2,sqrt((1-mu)*mu))
end

a = (1-mu)*mu*mu/sigma2-mu;
b = a*(1/mu-1);

%@test:1
%$ try
%$    [a, b] = beta_specification(.5, .05);
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    t(2) = abs(0.5-a/(a+b))<1e-12;
%$    t(3) = abs(0.05-a*b/((a+b)^2*(a+b+1)))<1e-12;
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ try
%$    [a, b] = beta_specification(0.5, .05, 1, 2);
%$    t(1) = false;
%$ catch
%$    t(1) = true;
%$ end
%$
%$ T = all(t);
%@eof:2

%@test:3
%$ try
%$    [a, b] = beta_specification(2.5, .05, 1, 2);
%$    t(1) = false;
%$ catch
%$    t(1) = true;
%$ end
%$
%$ T = all(t);
%@eof:3

%@test:4
%$ try
%$    [a, b] = beta_specification(.4, .6*.4+1e-4);
%$    t(1) = false;
%$ catch
%$    t(1) = true;
%$ end
%$
%$ T = all(t);
%@eof:4

%@test:5
%$ try
%$    [a, b] = beta_specification(.4, .6*.4+1e-6);
%$    t(1) = false;
%$ catch
%$    t(1) = true;
%$ end
%$ try
%$    [a, b] = beta_specification(.4, .6*.4-1e-6);
%$    t(2) = true;
%$ catch
%$    t(2) = false;
%$ end
%$
%$ T = all(t);
%@eof:5