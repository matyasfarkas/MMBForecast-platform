function [x,check] = solve1(func,x,j1,j2,jacobian_flag,gstep,tolf,tolx,maxit,debug,varargin)
% Solves systems of non linear equations of several variables
%
% INPUTS
%    func:            name of the function to be solved
%    x:               guess values
%    j1:              equations index for which the model is solved
%    j2:              unknown variables index
%    jacobian_flag=1: jacobian given by the 'func' function
%    jacobian_flag=0: jacobian obtained numerically
%    gstep            increment multiplier in numercial derivative
%                     computation
%    tolf             tolerance for residuals
%    tolx             tolerance for solution variation
%    maxit            maximum number of iterations
%    debug            debug flag
%    varargin:        list of extra arguments to the function
%
% OUTPUTS
%    x:               results
%    check=1:         the model can not be solved
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2001-2017 Dynare Team
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

nn = length(j1);
g = zeros(nn,1) ;

tolmin = tolx ;

stpmx = 100 ;

check = 0 ;

fvec = feval(func,x,varargin{:});
fvec = fvec(j1);

i = find(~isfinite(fvec));

if ~isempty(i)
    disp(['SOLVE1: during the resolution of the non-linear system, the evaluation of the following ' ...
          'equation(s) resulted in a non-finite number:'])
    disp(j1(i)')
    check = 1;
    return
end

f = 0.5*(fvec'*fvec) ;

if max(abs(fvec)) < tolf
    return ;
end

stpmax = stpmx*max([sqrt(x'*x);nn]) ;
first_time = 1;
if ~jacobian_flag
    fjac = zeros(nn,nn) ;
end
for its = 1:maxit
    if jacobian_flag
        [fvec,fjac] = feval(func,x,varargin{:});
        fvec = fvec(j1);
        fjac = fjac(j1,j2);
    else
        dh = max(abs(x(j2)),gstep(1)*ones(nn,1))*eps^(1/3);

        for j = 1:nn
            xdh = x ;
            xdh(j2(j)) = xdh(j2(j))+dh(j) ;
            t = feval(func,xdh,varargin{:});
            fjac(:,j) = (t(j1) - fvec)./dh(j) ;
            g(j) = fvec'*fjac(:,j) ;
        end
    end

    g = (fvec'*fjac)';
    if debug
        disp(['cond(fjac) ' num2str(condest(fjac))])
    end
    if issparse(fjac)
        rcond_fjac = 1/condest(fjac);
    else
        rcond_fjac = rcond(fjac);
    end
    if rcond_fjac < sqrt(eps)
        fjac2=fjac'*fjac;
        p=-(fjac2+sqrt(nn*eps)*max(sum(abs(fjac2)))*eye(nn))\(fjac'*fvec);
    else
        p = -fjac\fvec ;
    end
    xold = x ;
    fold = f ;

    [x, f, fvec, check] = lnsrch1(xold, fold, g, p, stpmax, func, j1, j2, tolx, varargin{:});

    if debug
        disp([its f])
        disp([xold x])
    end

    if check > 0
        den = max([f;0.5*nn]) ;
        if max(abs(g).*max([abs(x(j2)') ones(1,nn)])')/den < tolmin
            if max(abs(x(j2)-xold(j2))./max([abs(x(j2)') ones(1,nn)])') < tolx
                disp (' ')
                disp (['SOLVE: Iteration ' num2str(its)])
                disp (['Convergence on dX.'])
                disp (x)
                return
            end
        else
            disp (' ')
            disp (['SOLVE: Iteration ' num2str(its)])
            disp (['Spurious convergence.'])
            disp (x)
            return
        end

    elseif max(abs(fvec)) < tolf
        return
    end
end

check = 1;
skipline()
disp('SOLVE: maxit has been reached')

% 01/14/01 MJ lnsearch is now a separate function
% 01/16/01 MJ added varargin to function evaluation
% 04/13/01 MJ added test  f < tolf !!
% 05/11/01 MJ changed tests for 'check' so as to remove 'continue' which is
%             an instruction which appears only in version 6
