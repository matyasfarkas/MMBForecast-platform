function imf = zimpulse(Bh,swish,nn)
% Computing impulse functions with
%                imf = zimpulse(Bh,swish,nn)
%   where imf is in a format that is the SAME as in RATS.
%                Column: nvar responses to 1st shock,
%                            nvar responses to 2nd shock, and so on.
%                Row:  steps of impulse responses.
%         Bh is the estimated reduced form coefficient in the form
%              Y(T*nvar) = XB + U, X: T*k, B: k*nvar.  The matrix
%              form or dimension is the same as "Bh" from the function "sye";
%         swish is the inv(A0) in the structural model A(L)y(t) = e(t).
%         nn is the numbers of inputs [nvar,lags,# of impulse responses].
%  Written by Tao Zha
%
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


nvar = nn(1);
lags = nn(2);
imstep = nn(3);   % number of steps for impulse responses

Ah = Bh';
% Row: nvar equations
% Column: 1st lag (with nvar variables) to lags (with nvar variables) + const = k.

imf = zeros(imstep,nvar*nvar);
% Column: nvar responses to 1st shock, nvar responses to 2nd shock, and so on.
% Row:  steps of impulse responses.
M = zeros(nvar*(lags+1),nvar);
% Stack lags M's in the order of, e.g., [Mlags, ..., M2,M1;M0]
M(1:nvar,:) = swish;
Mtem = M(1:nvar,:);    % temporary M.
% first (initial) responses to 1 standard deviation shock.  Row: responses; Column: shocks
% * put in the form of "imf"
imf(1,:) = Mtem(:)';

t = 1;
ims1 = min([imstep-1 lags]);
while t <= ims1
   Mtem = Ah(:,1:nvar*t)*M(1:nvar*t,:);
   % Row: nvar equations, each for the nvar variables at tth lag
   M(nvar+1:nvar*(t+1),:)=M(1:nvar*t,:);
   M(1:nvar,:) = Mtem;
   imf(t+1,:) = Mtem(:)';
   % stack imf with each step, Row: 6 var to 1st shock, 6 var to 2nd shock, etc.
   t= t+1;
end

for t = lags+1:imstep-1
   Mtem = Ah(:,1:nvar*lags)*M(1:nvar*lags,:);
   % Row: nvar equations, each for the nvar variables at tth lag
   M(nvar+1:nvar*(t+1),:) = M(1:nvar*t,:);
   M(1:nvar,:)=Mtem;
   imf(t+1,:) = Mtem(:)';
   % stack imf with each step, Row: 6 var to 1st shock, 6 var to 2nd shock, etc.
end
