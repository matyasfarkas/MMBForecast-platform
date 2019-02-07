function H = bfgsi1(H0,dg,dx,Verbose,Save_files)
% H = bfgsi1(H0,dg,dx,Verbose,Save_files)
% Update Inverse Hessian matrix
%
% Inputs:
%   H0  [npar by npar]  initial inverse Hessian matrix
%   dg  [npar by 1]     previous change in gradient
%   dx  [npar by 1]     previous change in x;
%   Verbose [scalar]    Indicator for silent mode
%   Save_files [scalar] Indicator whether to save files
%
% 6/8/93 version that updates inverse Hessian instead of Hessian
% itself.
%
% Original file downloaded from:
% http://sims.princeton.edu/yftp/optimize/mfiles/bfgsi.m
%
% Copyright (C) 1993-2009 Christopher Sims
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

if size(dg,2)>1
    dg=dg';
end
if size(dx,2)>1
    dx=dx';
end
Hdg = H0*dg;
dgdx = dg'*dx;
if (abs(dgdx) >1e-12)
    H = H0 + (1+(dg'*Hdg)/dgdx)*(dx*dx')/dgdx - (dx*Hdg'+Hdg*dx')/dgdx;
else
    disp_verbose('bfgs update failed.',Verbose)
    disp_verbose(['|dg| = ' num2str(sqrt(dg'*dg)) '|dx| = ' num2str(sqrt(dx'*dx))],Verbose);
    disp_verbose(['dg''*dx = ' num2str(dgdx)],Verbose)
    disp_verbose(['|H*dg| = ' num2str(Hdg'*Hdg)],Verbose)
    H=H0;
end
if Save_files
    save('H.dat','H')
end
