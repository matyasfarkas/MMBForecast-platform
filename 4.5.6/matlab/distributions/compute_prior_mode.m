function m = compute_prior_mode(hyperparameters,shape) % --*-- Unitary tests --*--
% This function computes the mode of the prior distribution given the (two, three or four) hyperparameters
% of the prior distribution.
%
% INPUTS
%   hyperparameters     [double]    1*n vector of hyper parameters.
%   shape               [integer]   scalar specifying the prior shape:
%                                     shape=1 => Beta distribution,
%                                     shape=2 => Gamma distribution,
%                                     shape=3 => Gaussian distribution,
%                                     shape=4 => Inverse Gamma (type 1) distribution,
%                                     shape=5 => Uniform distribution,
%                                     shape=6 => Inverse Gamma (type 2) distribution,
%                                     shape=8 => Weibull distribution.
%
% OUTPUTS
%   m       [double]    scalar or 2*1 vector, the prior mode.
%
% REMARKS
% [1] The size of the vector of hyperparameters is 3 when the Gamma or Inverse Gamma is shifted and 4 when
%     the support of the Beta distribution is not [0,1].
% [2] The hyperparameters of the uniform distribution are the lower and upper bounds.
% [3] The uniform distribution has an infinity of modes. In this case the function returns the prior mean.
% [4] For the beta distribution we can have 1, 2 or an infinity of modes.

% Copyright (C) 2009-2017 Dynare Team
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

m = NaN ;
switch shape
  case 1
    if (hyperparameters(1)>1 && hyperparameters(2)>1)
        m = (hyperparameters(1)-1)/(hyperparameters(1)+hyperparameters(2)-2) ;
    elseif (hyperparameters(1)<1 && hyperparameters(2)<1)
        m = [ 0 ; 1 ] ;
    elseif ( hyperparameters(1)<1 && hyperparameters(2)>1-eps ) || ( abs(hyperparameters(1)-1)<2*eps && hyperparameters(2)>1 )
        m = 0;
    elseif ( hyperparameters(1)>1 && hyperparameters(2)<1+eps ) || ( abs(hyperparameters(1)-1)<2*eps && hyperparameters(2)<1 )
        m = 1;
    elseif ( abs(hyperparameters(1)-1)<2*eps && abs(hyperparameters(2)-1)<2*eps )% Uniform distribution!
        m = .5 ;
    end
    if length(hyperparameters)==4
        m = m*(hyperparameters(4)-hyperparameters(3)) + hyperparameters(3) ;
    end
  case 2
    % a = hyperparameters(1) [shape parameter]
    % b = hyperparameters(2) [scale parameter]
    if hyperparameters(1)<1
        m = 0;
    else
        m = (hyperparameters(1)-1)*hyperparameters(2);
    end
    if length(hyperparameters)>2
        m = m + hyperparameters(3);
    end
  case 3
    m = hyperparameters(1);
  case 4
    % s  = hyperparameters(1)
    % nu = hyperparameters(2)
    m = 1/sqrt((hyperparameters(2)+1)/hyperparameters(1));
    if length(hyperparameters)>2
        m = m + hyperparameters(3);
    end
  case 5
    m = hyperparameters(1);
  case 6
    % s  = hyperparameters(1)
    % nu = hyperparameters(2)
    m = hyperparameters(1)/(hyperparameters(2)+2) ;
    if length(hyperparameters)>2
        m = m + hyperparameters(3) ;
    end
  case 8
    % k = hyperparameters(1) [shape parameter]
    % s = hyperparameters(2) [scale parameter]
    if hyperparameters(1)<=1
        m = 0;
    else
        m = hyperparameters(2)*((hyperparameters(1)-1)/hyperparameters(1))^(1/hyperparameters(1));
    end
    if length(hyperparameters)>2
        % Add location parameter
        m = m + hyperparameters(3) ;
    end
  otherwise
    error('Unknown prior shape!')
end

%@test:1
%$ % Beta density
%$ try
%$     m1 = compute_prior_mode([2 1],1);
%$     m2 = compute_prior_mode([2 5 1 4],1); % Wolfram Alpha: BetaDistribution[2,5]
%$     t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ % Check the results
%$ if t(1)
%$     t(2) = dassert(m1,1,1e-6);
%$     t(3) = dassert(m2,1/5*3+1,1e-6);
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ % Gamma density
%$ try
%$     m1 = compute_prior_mode([1 2],2);
%$     m2 = compute_prior_mode([9 0.5 1],2);  % Wolfram Alpha: GammaDistribution[9,0.5]
%$     t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ % Check the results
%$ if t(1)
%$     t(2) = dassert(m1,0,1e-6);
%$     t(3) = dassert(m2,4+1,1e-6);
%$ end
%$ T = all(t);
%@eof:2

%@test:3
%$ % Normal density
%$ try
%$     m1 = compute_prior_mode([1 1],3);
%$     m2 = compute_prior_mode([2 5],3);
%$     t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ % Check the results
%$ if t(1)
%$     t(2) = dassert(m1,1,1e-6);
%$     t(3) = dassert(m2,2,1e-6);
%$ end
%$ T = all(t);
%@eof:3

%@test:4
%$ % Inverse Gamma I density
%$ try
%$     m1 = compute_prior_mode([8 2],4);
%$     m2 = compute_prior_mode([8 2 1],4);
%$     t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ % Check the results
%$ if t(1)
%$     t(2) = dassert(m1,1.632993161855452,1e-6);
%$     t(3) = dassert(m2,1.632993161855452+1,1e-6);
%$ end
%$ T = all(t);
%@eof:4

%@test:5
%$ % Uniform density
%$ try
%$     m1 = compute_prior_mode([0.5 1/sqrt(12)],5);
%$     m2 = compute_prior_mode([2 5 1 2],5);
%$     t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ % Check the results
%$ if t(1)
%$     t(2) = dassert(m1,0.5,1e-6);
%$     t(3) = dassert(m2,2,1e-6);
%$ end
%$ T = all(t);
%@eof:5

%@test:6
%$ % Inverse Gamma II density, parameterized with s and nu where  s=2*beta and nu=2*alpha
%$ try
%$     m1 = compute_prior_mode([8 2],6);  % Wolfram Alpha, parameterized with alpha beta: InversegammaDistribution[1,4]
%$     m2 = compute_prior_mode([8 4 1],6); % Wolfram Alpha, parameterized with alpha beta: InversegammaDistribution[2,4]
%$     t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ % Check the results
%$ if t(1)
%$     t(2) = dassert(m1,2,1e-6);
%$     t(3) = dassert(m2,1+4/3,1e-6);
%$ end
%$ T = all(t);
%@eof:6

%@test:7
%$ % Weibull density
%$ try
%$     m1 = compute_prior_mode([1 1],8);
%$     m2 = compute_prior_mode([2 1 1],8); % Wolfram Alpha: WeibullDistribution[2,1]
%$     t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ % Check the results
%$ if t(1)
%$     t(2) = dassert(m1,0,1e-6);
%$     t(3) = dassert(m2,1+1/sqrt(2),1e-6);
%$ end
%$ T = all(t);
%@eof:7

%@test:8
%$ % Unknown density
%$ try
%$     m1 = compute_prior_mode([1 1],7);
%$     t(1) = false;
%$ catch
%$     t(1) = true;
%$ end
%$
%$ T = all(t);
%@eof:8
