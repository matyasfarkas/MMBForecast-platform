function [logged_prior_density, dlprior, d2lprior, info] = priordens(x, pshape, p6, p7, p3, p4, initialization) % --*-- Unitary tests --*--
% Computes a prior density for the structural parameters of DSGE models
%
% INPUTS
%    x              [double]      vector with n elements.
%    pshape         [integer]     vector with n elements (bayestopt_.pshape).
%    p6:            [double]      vector with n elements, first  parameter of the prior distribution (bayestopt_.p6).
%    p7:            [double]      vector with n elements, second parameter of the prior distribution (bayestopt_.p7).
%    p3:            [double]      vector with n elements, lower bounds of the untruncated standard or generalized distribution
%    p4:            [double]      vector with n elements, upper bound of the untruncated standard or generalized distribution
%    initialization [integer]     if 1: initialize persistent variables
%
% OUTPUTS
%    logged_prior_density  [double]  scalar, log of the prior density evaluated at x.
%    info                  [double]  error code for index of Inf-prior parameter
%

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

persistent id1 id2 id3 id4 id5 id6 id8
persistent tt1 tt2 tt3 tt4 tt5 tt6 tt8

info=0;

if nargin > 6  && initialization
    % Beta indices.
    tt1 = true;
    id1 = find(pshape==1);
    if isempty(id1)
        tt1 = false;
    end
    % Gamma indices.
    tt2 = true;
    id2 = find(pshape==2);
    if isempty(id2)
        tt2 = false;
    end
    % Gaussian indices.
    tt3 = true;
    id3 = find(pshape==3);
    if isempty(id3)
        tt3 = false;
    end
    % Inverse-Gamma-1 indices.
    tt4 = true;
    id4 = find(pshape==4);
    if isempty(id4)
        tt4 = false;
    end
    % Uniform indices.
    tt5 = true;
    id5 = find(pshape==5);
    if isempty(id5)
        tt5 = false;
    end
    % Inverse-Gamma-2 indices.
    tt6 = true;
    id6 = find(pshape==6);
    if isempty(id6)
        tt6 = false;
    end
    % Weibull indices.
    tt8 = true;
    id8 = find(pshape==8);
    if isempty(id8)
        tt8 = false;
    end
end

logged_prior_density = 0.0;
dlprior = 0.0;
d2lprior = 0.0;

if tt1
    logged_prior_density = logged_prior_density + sum(lpdfgbeta(x(id1),p6(id1),p7(id1),p3(id1),p4(id1))) ;
    if isinf(logged_prior_density)
        if nargout ==4
            info=id1(isinf(lpdfgbeta(x(id1),p6(id1),p7(id1),p3(id1),p4(id1))));
        end
        return
    end
    if nargout == 2
        [tmp, dlprior(id1)]=lpdfgbeta(x(id1),p6(id1),p7(id1),p3(id1),p4(id1));
    elseif nargout == 3
        [tmp, dlprior(id1), d2lprior(id1)]=lpdfgbeta(x(id1),p6(id1),p7(id1),p3(id1),p4(id1));
    end
end

if tt2
    logged_prior_density = logged_prior_density + sum(lpdfgam(x(id2)-p3(id2),p6(id2),p7(id2))) ;
    if isinf(logged_prior_density)
        if nargout ==4
            info=id2(isinf(lpdfgam(x(id2)-p3(id2),p6(id2),p7(id2))));
        end
        return
    end
    if nargout == 2
        [tmp, dlprior(id2)]=lpdfgam(x(id2)-p3(id2),p6(id2),p7(id2));
    elseif nargout == 3
        [tmp, dlprior(id2), d2lprior(id2)]=lpdfgam(x(id2)-p3(id2),p6(id2),p7(id2));
    end
end

if tt3
    logged_prior_density = logged_prior_density + sum(lpdfnorm(x(id3),p6(id3),p7(id3))) ;
    if nargout == 2
        [tmp, dlprior(id3)]=lpdfnorm(x(id3),p6(id3),p7(id3));
    elseif nargout == 3
        [tmp, dlprior(id3), d2lprior(id3)]=lpdfnorm(x(id3),p6(id3),p7(id3));
    end
end

if tt4
    logged_prior_density = logged_prior_density + sum(lpdfig1(x(id4)-p3(id4),p6(id4),p7(id4))) ;
    if isinf(logged_prior_density)
        if nargout ==4
            info=id4(isinf(lpdfig1(x(id4)-p3(id4),p6(id4),p7(id4))));
        end
        return
    end
    if nargout == 2
        [tmp, dlprior(id4)]=lpdfig1(x(id4)-p3(id4),p6(id4),p7(id4));
    elseif nargout == 3
        [tmp, dlprior(id4), d2lprior(id4)]=lpdfig1(x(id4)-p3(id4),p6(id4),p7(id4));
    end
end

if tt5
    if any(x(id5)-p3(id5)<0) || any(x(id5)-p4(id5)>0)
        logged_prior_density = -Inf ;
        if nargout ==4
            info=id5((x(id5)-p3(id5)<0) || (x(id5)-p4(id5)>0));
        end
        return
    end
    logged_prior_density = logged_prior_density + sum(log(1./(p4(id5)-p3(id5)))) ;
    if nargout >1
        dlprior(id5)=zeros(length(id5),1);
    end
    if nargout == 3
        d2lprior(id5)=zeros(length(id5),1);
    end
end

if tt6
    logged_prior_density = logged_prior_density + sum(lpdfig2(x(id6)-p3(id6),p6(id6),p7(id6))) ;
    if isinf(logged_prior_density)
        if nargout ==4
            info=id6(isinf(lpdfig2(x(id6)-p3(id6),p6(id6),p7(id6))));
        end
        return
    end
    if nargout == 2
        [tmp, dlprior(id6)]=lpdfig2(x(id6)-p3(id6),p6(id6),p7(id6));
    elseif nargout == 3
        [tmp, dlprior(id6), d2lprior(id6)]=lpdfig2(x(id6)-p3(id6),p6(id6),p7(id6));
    end
end

if tt8
    logged_prior_density = logged_prior_density + sum(lpdfgweibull(x(id8),p6(id8),p7(id8)));
    if isinf(logged_prior_density)
        if nargout ==4
            info=id8(isinf(log(lpdfgweibull(x(id8),p6(id8),p7(id8)))));
        end
        return
    end
    if nargout==2
        [tmp, dlprior(id8)] = lpdfgweibull(x(id8),p6(id8),p7(id8))
    elseif nargout==3
        [tmp, dlprior(id8), ds2lprior(id8)] = lpdfgweibull(x(id8),p6(id8),p7(id8))
    end
end

if nargout==3
    d2lprior = diag(d2lprior);
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
%$ % Call the tested routine
%$ try
%$     % Initialization of priordens.
%$     lpdstar = priordens(p5, p0, p6, p7, p3, p4, true);
%$     % Do simulations in a loop and estimate recursively the mean and teh variance.
%$     LPD = NaN(10000,1);
%$     for i = 1:10000
%$          draw = p5+randn(size(p5))*.02;
%$          LPD(i) = priordens(p5, p0, p6, p7, p3, p4);
%$     end
%$     t(1) = true;
%$ catch
%$     t(1) = false;
%$ end
%$
%$ if t(1)
%$     t(2) = all(LPD<=lpdstar);
%$ end
%$ T = all(t);
%@eof:1