function hessian_mat = hessian_sparse(func,x,gstep,varargin)
% function hessian_mat = hessian_sparse(func,x,gstep,varargin)
% Computes second order partial derivatives
%
% INPUTS
%    func        [string]   name of the function
%    x           [double]   vector, the Hessian of "func" is evaluated at x.
%    gstep       [double]   scalar, size of epsilon.
%    varargin    [void]     list of additional arguments for "func".
%
% OUTPUTS
%    hessian_mat [double, sparse]   Hessian matrix
%
% ALGORITHM
%    Uses Abramowitz and Stegun (1965) formulas 25.3.24 and 25.3.27 p. 884
%
% SPECIAL REQUIREMENTS
%    none
%

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

if ~isa(func, 'function_handle')
    func = str2func(func);
end
n=size(x,1);
h1=max(abs(x),sqrt(gstep(1))*ones(n,1))*eps^(1/6)*gstep(2);
h_1=h1;
xh1=x+h1;
h1=xh1-x;
xh1=x-h_1;
h_1=x-xh1;
xh1=x;
f0=feval(func,x,varargin{:});
f1=zeros(size(f0,1),n);
f_1=f1;
for i=1:n
    xh1(i)=x(i)+h1(i);
    f1(:,i)=feval(func,xh1,varargin{:});
    xh1(i)=x(i)-h_1(i);
    f_1(:,i)=feval(func,xh1,varargin{:});
    xh1(i)=x(i);
end
xh_1=xh1;
hessian_mat = sparse(size(f0,1),n*n);

for i=1:n
    %     if i > 1
    %         k=[i:n:n*(i-1)];
    %         hessian_mat(:,(i-1)*n+1:(i-1)*n+i-1)=hessian_mat(:,k);
    %         hessian_mat(:,k)=0;
    %     end
    hessian_mat(:,(i-1)*n+i)=(f1(:,i)+f_1(:,i)-2*f0)./(h1(i)*h_1(i));
    temp=f1+f_1-f0*ones(1,n);
    for j=1:i-1
        xh1(i)=x(i)+h1(i);
        xh1(j)=x(j)+h_1(j);
        xh_1(i)=x(i)-h1(i);
        xh_1(j)=x(j)-h_1(j);
        hessian_mat(:,(i-1)*n+j)=-(-feval(func,xh1,varargin{:})-feval(func,xh_1,varargin{:})+temp(:,i)+temp(:,j))./(2*h1(i)*h_1(j));
        xh1(i)=x(i);
        xh1(j)=x(j);
        xh_1(i)=x(i);
        xh_1(j)=x(j);
    end
end