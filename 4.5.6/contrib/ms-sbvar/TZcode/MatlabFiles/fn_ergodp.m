function gpi = fn_ergodp(P)
% gpi = fn_ergodp(P)
%    Compute the ergodic probabilities.  See Hamilton p.681.
%
% P:  n-by-n matrix of transition matrix where all elements in each column sum up to 1.
%-----
% gpi:  n-by-1 vector of ergodic probabilities.
%
% Tao Zha August 2000
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


[gpim,gpid] = eig(P);  % m: matrix; d: diagonal
[gpidv,gpidvinx] = sort(diag(gpid));
gpidv = fliplr(gpidv);
gpidvinx = flipud(gpidvinx);
gpim = gpim(:,gpidvinx);
gpi = gpim(:,1)/sum(gpim(:,1));
