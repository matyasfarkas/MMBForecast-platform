function xcum = priorcdf(para, pshape, p6, p7, p3, p4)

% This procedure transforms x vectors into cumulative values
% pshape: 0 is point mass, both para and p2 are ignored
%         1 is BETA(mean,stdd)
%         2 is GAMMA(mean,stdd)
%         3 is NORMAL(mean,stdd)
%         4 is INVGAMMA(s^2,nu) type I
%         5 is UNIFORM [p1,p2]
%         6 is INNGAMMA(s^2,nu) type II
%         8 is WEIBULL(s, k)
% Adapted by M. Ratto from MJ priordens.m

% Copyright (C) 2012-2015 Dynare Team
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

for i=1:length(pshape)
    switch pshape(i)
      case 1 % (generalized) BETA Prior
        para(:,i) = (para(:,i)-p3(i))./(p4(i)-p3(i));
        xcum(:,i) = betainc(para(:,i),p6(i),p7(i));
      case 2 % GAMMA PRIOR
        xcum(:,i) = gamcdf(para(:,i)-p3(i),p6(i),p7(i));
      case 3 % GAUSSIAN PRIOR
        xcum(:,i) = 0.5 * erfc(-(para(:,i)-p6(i))/p7(i) ./ sqrt(2));
      case 4 % INVGAMMA1 PRIOR
        xcum(:,i) = gamcdf(1./(para(:,i)-p3(i)).^2,p7(i)/2,2/p6(i));
      case 5 % UNIFORM PRIOR
        xcum(:,i) = (para(:,i)-p3(i))./(p4(i)-p3(i));
      case 6 % INVGAMMA2 PRIOR
        xcum(:,i) = gamcdf(1./(para(:,i)-p3(i)),p7(i)/2,2/p6(i));
      case 8 % WEIBULL
        xcum(:,i) = wblcdf(para(:,i)-p3(i),p6(i),p7(i));
      otherwise
        error('Unknown prior shape!')
    end
end
