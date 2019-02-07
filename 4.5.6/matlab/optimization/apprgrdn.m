function g = apprgrdn(x,f,fun,deltax,obj,varargin)
% g = apprgrdn(x,f,fun,deltax,obj,varargin)
% Performs the finite difference approximation of the gradient <g> at a
% point <x> used in solveopt
%
% Inputs:
% x:        point at which to evaluate gradient
% f:        calculated function value at a point x;
% fun:      Name of the Matlab function calculating the function values
% deltax:   vector of the relative stepsizes,
% obj       flag indicating whether the gradient of the objective
%           function (1) or the constraint function (0) is to be calculated.
%
% Modified by Giovanni Lombardo and Johannes Pfeifer to accomodate Dynare
% structure
%
%
% Copyright (C) 1997-2008, Alexei Kuntsevich and Franz Kappel
% Copyright (C) 2008-2015 Giovanni Lombardo
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

n=max(size(x)); ee=ones(size(x));
di=abs(x); idx=find(di<5e-15); di(idx)=5e-15*ee(idx);
di=deltax.*di;
if obj
    idx=find(abs(di)<2e-10);
    di(idx)=2e-10*sign(di(idx));
else
    idx=find(abs(di)<5e-15);
    di(idx)=5e-15*sign(di(idx));
end
y=x;

g=NaN(n,1);
for i=1:n
    y(i)=x(i)+di(i);
    fi=feval(fun,y,varargin{:});
    if obj
        if fi==f
            for j=1:3
                di(i)=di(i)*10;  y(i)=x(i)+di(i);
                fi=feval(fun,y,varargin{:});
                if fi~=f
                    break
                end
            end
        end
    end
    g(i)=(fi-f)/di(i);
    if obj
        if ~isempty(idx) && any(idx==i)
            y(i)=x(i)-di(i);
            fi=feval(fun,y,varargin{:});
            g(i)=.5*(g(i)+(f-fi)/di(i));
        end
    end
    y(i)=x(i);
end
