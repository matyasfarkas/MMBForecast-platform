function [shape, scale] = weibull_specification(mu, sigma2, lb, name)   % --*-- Unitary tests --*--

% Returns the hyperparameters of the Weibull distribution given the expectation and variance.
%
% INPUTS
%
%
% OUTPUTS
%
%

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
    error('The prior expectation (%f) %scannot be smaller than the lower bound of the Weibull distribution (%f)!', mu, name1, lb)
end

if isinf(sigma2)
    error('The variance of the Gamma distribution has to be finite%s!', name2)
end

scale = NaN;
shape = NaN;

mu = mu-lb;
mu2 = mu*mu;

eqn = @(k) gammaln(1+2./k) - 2*gammaln(1+1./k) - log(1+sigma2/mu2);
eqn2 = @(k) eqn(k).*eqn(k);

kstar = fminbnd(eqn2, 1e-9, 100);
[shape, fval, exitflag]  = fzero(eqn, kstar);

if exitflag<1
    shape = NaN;
    return
end

scale = mu/gamma(1+1/shape);

%@test:1
%$ debug = false;
%$ scales = 1:.01:5;
%$ shapes = .5:.01:2;
%$ n_scales = length(scales);
%$ n_shapes = length(shapes);
%$ mu = NaN(n_scales, n_shapes);
%$ s2 = NaN(n_scales, n_shapes);
%$ for i=1:n_shapes
%$     g1 = gamma(1+1/shapes(i));
%$     g2 = gamma(1+2/shapes(i));
%$     g3 = g1*g1;
%$     for j=1:n_scales
%$         mu(j, i) = scales(j)*g1;
%$         s2(j, i) = scales(j)*scales(j)*(g2-g3);
%$     end
%$ end
%$ if debug
%$    success = [];
%$    failed1 = [];
%$    failed1_ = [];
%$    failed2 = [];
%$ end
%$ try
%$    for i=1:n_shapes
%$        for j=1:n_scales
%$           if debug
%$               disp(sprintf('... mu=%s and s2=%s', num2str(mu(j,i)),num2str(s2(j,i))))
%$           end
%$           if ~isnan(mu(j,i)) && ~isnan(s2(j,i)) && ~isinf(mu(j,i)) && ~isinf(s2(j,i))
%$               [shape, scale] = weibull_specification(mu(j,i), s2(j,i));
%$               if isnan(scale)
%$                  t = false;
%$               else
%$                  if abs(scales(j)-scale)<1e-9 && abs(shapes(i)-shape)<1e-9
%$                      t = true;
%$                  else
%$                      t = false;
%$                  end
%$               end
%$               if ~t && debug
%$                  failed1 = [failed1; mu(j,i) s2(j,i)];
%$                  failed1_ = [failed1_; shapes(i) scales(j)];
%$                  error('UnitTest','Cannot compute scale and shape hyperparameters for mu=%s and s2=%s', num2str(mu(j,i)), num2str(s2(j,i)))
%$               end
%$               if debug
%$                   success = [success; mu(j,i) s2(j,i)];
%$               end
%$           else
%$               failed2 = [failed2; shapes(i) scales(j)];
%$               continue % Pass this test
%$           end
%$       end
%$   end
%$ catch
%$    t = false;
%$ end
%$
%$ if debug
%$     figure(1)
%$     plot(success(:,1),success(:,2),'ok');
%$     if ~isempty(failed1)
%$         hold on
%$         plot(failed1(:,1),failed1(:,2),'or');
%$         hold off
%$         figure(2)
%$         plot(failed1_(:,1),failed1_(:,2),'or')
%$     end
%$     if ~isempty(failed2)
%$         figure(2)
%$         plot(failed2(:,1),failed2(:,2),'or');
%$     end
%$ end
%$ T = all(t);
%@eof:1
