function pdraw = prior_draw(BayesInfo, prior_trunc, uniform) % --*-- Unitary tests --*--

% This function generate one draw from the joint prior distribution and
% allows sampling uniformly from the prior support (uniform==1 when called with init==1)
%
% INPUTS
%   o init             [integer]    scalar equal to:
%                                       1: first call to set up persistent variables
%                                             describing the prior
%                                       0: subsequent call to get prior
%                                               draw
%   o uniform          [integer]    scalar used in initialization (init=1), equal to:
%                                       1: sample uniformly from prior
%                                           support (overwrites prior shape used for sampling within this function)
%                                       0: sample from joint prior distribution
%
% OUTPUTS
%   o pdraw            [double]     1*npar vector, draws from the joint prior density.
%
%
% SPECIAL REQUIREMENTS
%   none
%
% NOTE 1. Input arguments 1 an 2 are only needed for initialization.
% NOTE 2. A given draw from the joint prior distribution does not satisfy BK conditions a priori.
% NOTE 3. This code relies on bayestopt_ as created in the base workspace
%           by the preprocessor (or as updated in subsequent pieces of code and handed to the base workspace)
%
% Copyright (C) 2006-2017 Dynare Team
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

persistent p6 p7 p3 p4 lb ub
persistent uniform_index gaussian_index gamma_index beta_index inverse_gamma_1_index inverse_gamma_2_index weibull_index
persistent uniform_draws gaussian_draws gamma_draws beta_draws inverse_gamma_1_draws inverse_gamma_2_draws weibull_draws

if nargin>0
    p6 = BayesInfo.p6;
    p7 = BayesInfo.p7;
    p3 = BayesInfo.p3;
    p4 = BayesInfo.p4;
    bounds = prior_bounds(BayesInfo, prior_trunc);
    lb = bounds.lb;
    ub = bounds.ub;
    number_of_estimated_parameters = length(p6);
    if nargin>2 && uniform
        prior_shape = repmat(5,number_of_estimated_parameters,1);
    else
        prior_shape = BayesInfo.pshape;
    end
    beta_index = find(prior_shape==1);
    if isempty(beta_index)
        beta_draws = false;
    else
        beta_draws = true;
    end
    gamma_index = find(prior_shape==2);
    if isempty(gamma_index)
        gamma_draws = false;
    else
        gamma_draws = true;
    end
    gaussian_index = find(prior_shape==3);
    if isempty(gaussian_index)
        gaussian_draws = false;
    else
        gaussian_draws = true;
    end
    inverse_gamma_1_index = find(prior_shape==4);
    if isempty(inverse_gamma_1_index)
        inverse_gamma_1_draws = false;
    else
        inverse_gamma_1_draws = true;
    end
    uniform_index = find(prior_shape==5);
    if isempty(uniform_index)
        uniform_draws = false;
    else
        uniform_draws = true;
    end
    inverse_gamma_2_index = find(prior_shape==6);
    if isempty(inverse_gamma_2_index)
        inverse_gamma_2_draws = false;
    else
        inverse_gamma_2_draws = true;
    end
    weibull_index = find(prior_shape==8);
    if isempty(weibull_index)
        weibull_draws = false;
    else
        weibull_draws = true;
    end
    pdraw = NaN(number_of_estimated_parameters,1);
    return
end

if uniform_draws
    pdraw(uniform_index) = rand(length(uniform_index),1).*(p4(uniform_index)-p3(uniform_index)) + p3(uniform_index);
    out_of_bound = find( (pdraw(uniform_index)'>ub(uniform_index)) | (pdraw(uniform_index)'<lb(uniform_index)));
    while ~isempty(out_of_bound)
        pdraw(uniform_index) = rand(length(uniform_index),1).*(p4(uniform_index)-p3(uniform_index)) + p3(uniform_index);
        out_of_bound = find( (pdraw(uniform_index)'>ub(uniform_index)) | (pdraw(uniform_index)'<lb(uniform_index)));
    end
end

if gaussian_draws
    pdraw(gaussian_index) = randn(length(gaussian_index),1).*p7(gaussian_index) + p6(gaussian_index);
    out_of_bound = find( (pdraw(gaussian_index)'>ub(gaussian_index)) | (pdraw(gaussian_index)'<lb(gaussian_index)));
    while ~isempty(out_of_bound)
        pdraw(gaussian_index(out_of_bound)) = randn(length(gaussian_index(out_of_bound)),1).*p7(gaussian_index(out_of_bound)) + p6(gaussian_index(out_of_bound));
        out_of_bound = find( (pdraw(gaussian_index)'>ub(gaussian_index)) | (pdraw(gaussian_index)'<lb(gaussian_index)));
    end
end

if gamma_draws
    pdraw(gamma_index) = gamrnd(p6(gamma_index),p7(gamma_index))+p3(gamma_index);
    out_of_bound = find( (pdraw(gamma_index)'>ub(gamma_index)) | (pdraw(gamma_index)'<lb(gamma_index)));
    while ~isempty(out_of_bound)
        pdraw(gamma_index(out_of_bound)) = gamrnd(p6(gamma_index(out_of_bound)),p7(gamma_index(out_of_bound)))+p3(gamma_index(out_of_bound));
        out_of_bound = find( (pdraw(gamma_index)'>ub(gamma_index)) | (pdraw(gamma_index)'<lb(gamma_index)));
    end
end

if beta_draws
    pdraw(beta_index) = (p4(beta_index)-p3(beta_index)).*betarnd(p6(beta_index),p7(beta_index))+p3(beta_index);
    out_of_bound = find( (pdraw(beta_index)'>ub(beta_index)) | (pdraw(beta_index)'<lb(beta_index)));
    while ~isempty(out_of_bound)
        pdraw(beta_index(out_of_bound)) = (p4(beta_index(out_of_bound))-p3(beta_index(out_of_bound))).*betarnd(p6(beta_index(out_of_bound)),p7(beta_index(out_of_bound)))+p3(beta_index(out_of_bound));
        out_of_bound = find( (pdraw(beta_index)'>ub(beta_index)) | (pdraw(beta_index)'<lb(beta_index)));
    end
end

if inverse_gamma_1_draws
    pdraw(inverse_gamma_1_index) = ...
        sqrt(1./gamrnd(p7(inverse_gamma_1_index)/2,2./p6(inverse_gamma_1_index)))+p3(inverse_gamma_1_index);
    out_of_bound = find( (pdraw(inverse_gamma_1_index)'>ub(inverse_gamma_1_index)) | (pdraw(inverse_gamma_1_index)'<lb(inverse_gamma_1_index)));
    while ~isempty(out_of_bound)
        pdraw(inverse_gamma_1_index(out_of_bound)) = ...
            sqrt(1./gamrnd(p7(inverse_gamma_1_index(out_of_bound))/2,2./p6(inverse_gamma_1_index(out_of_bound))))+p3(inverse_gamma_1_index(out_of_bound));
        out_of_bound = find( (pdraw(inverse_gamma_1_index)'>ub(inverse_gamma_1_index)) | (pdraw(inverse_gamma_1_index)'<lb(inverse_gamma_1_index)));
    end
end

if inverse_gamma_2_draws
    pdraw(inverse_gamma_2_index) = ...
        1./gamrnd(p7(inverse_gamma_2_index)/2,2./p6(inverse_gamma_2_index))+p3(inverse_gamma_2_index);
    out_of_bound = find( (pdraw(inverse_gamma_2_index)'>ub(inverse_gamma_2_index)) | (pdraw(inverse_gamma_2_index)'<lb(inverse_gamma_2_index)));
    while ~isempty(out_of_bound)
        pdraw(inverse_gamma_2_index(out_of_bound)) = ...
            1./gamrnd(p7(inverse_gamma_2_index(out_of_bound))/2,2./p6(inverse_gamma_2_index(out_of_bound)))+p3(inverse_gamma_2_index(out_of_bound));
        out_of_bound = find( (pdraw(inverse_gamma_2_index)'>ub(inverse_gamma_2_index)) | (pdraw(inverse_gamma_2_index)'<lb(inverse_gamma_2_index)));
    end
end

if weibull_draws
    pdraw(weibull_index) = wblrnd(p7(weibull_index), p6(weibull_index)) + p3(weibull_index);
    out_of_bound = find( (pdraw(weibull_index)'>ub(weibull_index)) | (pdraw(weibull_index)'<lb(weibull_index)));
    while ~isempty(out_of_bound)
        pdraw(weibull_index(out_of_bound)) = wblrnd(p7(weibull_index(out_of_bound)),p6(weibull_index(out_of_bound)))+p3(weibull_index(out_of_bound));
        out_of_bound = find( (pdraw(weibull_index)'>ub(weibull_index)) | (pdraw(weibull_index)'<lb(weibull_index)));
    end
end

%@test:1
%$ % Fill global structures with required fields...
%$ prior_trunc = 1e-10;
%$ p0 = repmat([1; 2; 3; 4; 5; 6; 8], 2, 1);    % Prior shape
%$ p1 = .4*ones(14,1);                          % Prior mean
%$ p2 = .2*ones(14,1);                          % Prior std.
%$ p3 = NaN(14,1);
%$ p4 = NaN(14,1);
%$ p5 = NaN(14,1);
%$ p6 = NaN(14,1);
%$ p7 = NaN(14,1);
%$
%$ for i=1:14
%$     switch p0(i)
%$       case 1
%$         % Beta distribution
%$         p3(i) = 0;
%$         p4(i) = 1;
%$         [p6(i), p7(i)] = beta_specification(p1(i), p2(i)^2, p3(i), p4(i));
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 1);
%$       case 2
%$         % Gamma distribution
%$         p3(i) = 0;
%$         p4(i) = Inf;
%$         [p6(i), p7(i)] = gamma_specification(p1(i), p2(i)^2, p3(i), p4(i));
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 2);
%$       case 3
%$         % Normal distribution
%$         p3(i) = -Inf;
%$         p4(i) = Inf;
%$         p6(i) = p1(i);
%$         p7(i) = p2(i);
%$         p5(i) = p1(i);
%$       case 4
%$         % Inverse Gamma (type I) distribution
%$         p3(i) = 0;
%$         p4(i) = Inf;
%$         [p6(i), p7(i)] = inverse_gamma_specification(p1(i), p2(i)^2, p3(i), 1, false);
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 4);
%$       case 5
%$         % Uniform distribution
%$         [p1(i), p2(i), p6(i), p7(i)] = uniform_specification(p1(i), p2(i), p3(i), p4(i));
%$         p3(i) = p6(i);
%$         p4(i) = p7(i);
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 5);
%$       case 6
%$         % Inverse Gamma (type II) distribution
%$         p3(i) = 0;
%$         p4(i) = Inf;
%$         [p6(i), p7(i)] = inverse_gamma_specification(p1(i), p2(i)^2, p3(i), 2, false);
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 6);
%$       case 8
%$         % Weibull distribution
%$         p3(i) = 0;
%$         p4(i) = Inf;
%$         [p6(i), p7(i)] = weibull_specification(p1(i), p2(i)^2, p3(i));
%$         p5(i) = compute_prior_mode([p6(i) p7(i)], 8);
%$       otherwise
%$         error('This density is not implemented!')
%$     end
%$ end
%$
%$ BayesInfo.pshape = p0;
%$ BayesInfo.p1 = p1;
%$ BayesInfo.p2 = p2;
%$ BayesInfo.p3 = p3;
%$ BayesInfo.p4 = p4;
%$ BayesInfo.p5 = p5;
%$ BayesInfo.p6 = p6;
%$ BayesInfo.p7 = p7;
%$
%$ ndraws = 1e5;
%$ m0 = BayesInfo.p1; %zeros(14,1);
%$ v0 = diag(BayesInfo.p2.^2); %zeros(14);
%$
%$ % Call the tested routine
%$ try
%$     % Initialization of prior_draws.
%$     prior_draw(BayesInfo, prior_trunc, false);
%$     % Do simulations in a loop and estimate recursively the mean and teh variance.
%$     for i = 1:ndraws
%$          draw = transpose(prior_draw());
%$          m1 = m0 + (draw-m0)/i;
%$          m2 = m1*m1';
%$          v0 = v0 + ((draw*draw'-m2-v0) + (i-1)*(m0*m0'-m2'))/i;
%$          m0 = m1;
%$     end
%$     t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ if t(1)
%$     t(2) = all(abs(m0-BayesInfo.p1)<3e-3);
%$     t(3) = all(all(abs(v0-diag(BayesInfo.p2.^2))<2e-3));
%$ end
%$ T = all(t);
%@eof:1