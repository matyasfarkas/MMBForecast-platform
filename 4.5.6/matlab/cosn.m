function [co, b, yhat] = cosn(H)

% function co = cosn(H);
% computes the cosine of the angle between the H(:,1) and its
% projection onto the span of H(:,2:end)
%
% Not the same as multiple correlation coefficient since the means are not
% zero
%
% Copyright (C) 2008-2017 Dynare Team
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

y = H(:,1);
X = H(:,2:end);

b=(X\y);
if any(isnan(b)) || any(isinf(b))
    b=0;
end
yhat =  X*b;
if rank(yhat)
    co = abs(y'*yhat/sqrt((y'*y)*(yhat'*yhat)));
else
    co=0;
end
