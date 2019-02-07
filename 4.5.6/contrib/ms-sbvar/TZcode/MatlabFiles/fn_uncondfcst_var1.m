function y = fn_uncondfcst_var1(G1, y0, nsteps);
%y = fn_irf_var1(G1,y0,nsteps);
% Inputs:
%  G1: n-by-n;
%  y0: n-by-1, initial condition;
%  nsteps: number of forecasts steps.
%---
% Outputs:
%  y: nsteps-by-n unconditional forecasts.
%
% See fn_vds.m, fn_irf_var1.m.
%
% Copyright (C) 1997-2012 Tao Zha
%
% This free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% It is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% If you did not received a copy of the GNU General Public License
% with this software, see <http://www.gnu.org/licenses/>.
%

n = length(y0);

[n1,n2] = size(G1);
if (n1 ~= n2) || (n1 ~= n)
   error('fn_uncondfcst_var1.m: make sure that (1) G1 is square and (2) size(G1,1) = length(y0)');
end

y = zeros(nsteps,n);


%---- Forecast at the first step.
y(1,:) = (G1*y0)';

for ti = 2:nsteps
   y(ti,:) = (G1*y(ti-1,:)')';
end
