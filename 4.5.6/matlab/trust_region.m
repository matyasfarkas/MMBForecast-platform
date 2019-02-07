function [x,check,info] = trust_region(fcn,x0,j1,j2,jacobian_flag,gstep,tolf,tolx,maxiter,debug,varargin)
% Solves systems of non linear equations of several variables, using a
% trust-region method.
%
% INPUTS
%    fcn:             name of the function to be solved
%    x0:              guess values
%    j1:              equations index for which the model is solved
%    j2:              unknown variables index
%    jacobian_flag=1: jacobian given by the 'func' function
%    jacobian_flag=0: jacobian obtained numerically
%    gstep            increment multiplier in numercial derivative
%                     computation
%    tolf             tolerance for residuals
%    tolx             tolerance for solution variation
%    maxiter          maximum number of iterations
%    debug            debug flag
%    varargin:        list of arguments following bad_cond_flag
%
% OUTPUTS
%    x:               results
%    check=1:         the model can not be solved
%    info:            detailed exitcode
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2008-2012 VZLU Prague, a.s.
% Copyright (C) 2014-2017 Dynare Team
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
%
% Initial author: Jaroslav Hajek <highegg@gmail.com>, for GNU Octave

if (ischar (fcn))
    fcn = str2func (fcn);
end

n = length(j1);

% These defaults are rather stringent. I think that normally, user
% prefers accuracy to performance.

macheps = eps (class (x0));

niter = 1;

x = x0;
info = 0;

% Initial evaluation.
% Handle arbitrary shapes of x and f and remember them.
fvec = fcn (x, varargin{:});
fvec = fvec(j1);
fn = norm (fvec);

% Outer loop.
while (niter < maxiter && ~info)

    % Calculate function value and Jacobian (possibly via FD).
    if jacobian_flag
        [fvec, fjac] = fcn (x, varargin{:});
        fvec = fvec(j1);
        fjac = fjac(j1,j2);
    else
        dh = max(abs(x(j2)),gstep(1)*ones(n,1))*eps^(1/3);

        for j = 1:n
            xdh = x ;
            xdh(j2(j)) = xdh(j2(j))+dh(j) ;
            t = fcn(xdh,varargin{:});
            fjac(:,j) = (t(j1) - fvec)./dh(j) ;
        end
    end

    % Get column norms, use them as scaling factors.
    jcn = sqrt(sum(fjac.*fjac))';
    if (niter == 1)
        dg = jcn;
        dg(dg == 0) = 1;
    else
        % Rescale adaptively.
        % FIXME: the original minpack used the following rescaling strategy:
        %   dg = max (dg, jcn);
        % but it seems not good if we start with a bad guess yielding Jacobian
        % columns with large norms that later decrease, because the corresponding
        % variable will still be overscaled. So instead, we only give the old
        % scaling a small momentum, but do not honor it.

        dg = max (0.1*dg, jcn);
    end

    if (niter == 1)
        xn = norm (dg .* x(j2));
        % FIXME: something better?
        delta = max (xn, 1);
    end

    % Get trust-region model (dogleg) minimizer.
    s = - dogleg (fjac, fvec, dg, delta);
    w = fvec + fjac * s;

    sn = norm (dg .* s);
    if (niter == 1)
        delta = min (delta, sn);
    end

    x2 = x;
    x2(j2) = x2(j2) + s;
    fvec1 = fcn (x2, varargin{:});
    fvec1 = fvec1(j1);
    fn1 = norm (fvec1);

    if (fn1 < fn)
        % Scaled actual reduction.
        actred = 1 - (fn1/fn)^2;
    else
        actred = -1;
    end

    % Scaled predicted reduction, and ratio.
    t = norm (w);
    if (t < fn)
        prered = 1 - (t/fn)^2;
        ratio = actred / prered;
    else
        prered = 0;
        ratio = 0;
    end

    % Update delta.
    if (ratio < 0.1)
        delta = 0.5*delta;
        if (delta <= 1e1*macheps*xn)
            % Trust region became uselessly small.
            if (fn1 <= tolf)
                info = 1;
            else
                info = -3;
            end
            break
        end
    elseif (abs (1-ratio) <= 0.1)
        delta = 1.4142*sn;
    elseif (ratio >= 0.5)
        delta = max (delta, 1.4142*sn);
    end

    if (ratio >= 1e-4)
        % Successful iteration.
        x(j2) = x(j2) + s;
        xn = norm (dg .* x(j2));
        fvec = fvec1;
        fn = fn1;
    end

    niter = niter + 1;


    % Tests for termination conditions. A mysterious place, anything
    % can happen if you change something here...

    % The rule of thumb (which I'm not sure M*b is quite following)
    % is that for a tolerance that depends on scaling, only 0 makes
    % sense as a default value. But 0 usually means uselessly long
    % iterations, so we need scaling-independent tolerances wherever
    % possible.

    if (fn <= tolf)
        info = 1;
    end
end
if info==1
    check = 0;
else
    check = 1;
end
end


% Solve the double dogleg trust-region least-squares problem:
% Minimize norm(r*x-b) subject to the constraint norm(d.*x) <= delta,
% x being a convex combination of the gauss-newton and scaled gradient.

% TODO: error checks
% TODO: handle singularity, or leave it up to mldivide?

function x = dogleg (r, b, d, delta)
% Get Gauss-Newton direction.
x = r \ b;
xn = norm (d .* x);
if (xn > delta)
    % GN is too big, get scaled gradient.
    s = (r' * b) ./ d;
    sn = norm (s);
    if (sn > 0)
        % Normalize and rescale.
        s = (s / sn) ./ d;
        % Get the line minimizer in s direction.
        tn = norm (r*s);
        snm = (sn / tn) / tn;
        if (snm < delta)
            % Get the dogleg path minimizer.
            bn = norm (b);
            dxn = delta/xn; snmd = snm/delta;
            t = (bn/sn) * (bn/xn) * snmd;
            t = t - dxn * snmd^2 + sqrt ((t-dxn)^2 + (1-dxn^2)*(1-snmd^2));
            alpha = dxn*(1-snmd^2) / t;
        else
            alpha = 0;
        end
    else
        alpha = delta / xn;
        snm = 0;
    end
    % Form the appropriate convex combination.
    x = alpha * x + ((1-alpha) * min (snm, delta)) * s;
end
end
