function [Da1,DP1,D2a,D2P] = univariate_computeDstate(k,a,P,T,Da,DP,DT,DOm,notsteady,D2a,D2P,D2T,D2Om)

% Copyright (C) 2012-2017 Dynare Team
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

% AUTHOR(S) marco.ratto@jrc.ec.europa.eu


DP1=DP*0;
Da1=Da*0;
for j=1:k
    Da1(:,j) = T*Da(:,j) + DT(:,:,j)*a;
    if notsteady
        DP1(:,:,j) = T*DP(:,:,j)*T'+DT(:,:,j)*P*T'+T*P*DT(:,:,j)';
    else
        DP1=DP;
    end
end
if notsteady
    DP1 = DP1 + DOm;
end
if nargout>2
    jcount=0;
    for j=1:k
        for i=1:j
            jcount=jcount+1;
            D2a(:,j,i) = DT(:,:,i)*Da(:,j) + DT(:,:,j)*Da(:,i) + T*D2a(:,j,i)+ reshape(D2T(:,jcount),size(T))*a;
            D2a(:,i,j) = D2a(:,j,i);
            if notsteady
                tmp = dyn_unvech(D2P(:,jcount));
                tmp = T*tmp*T' +DT(:,:,i)*DP(:,:,j)*T'+T*DP(:,:,j)*DT(:,:,i)' + ...
                      DT(:,:,j)*DP(:,:,i)*T'+T*DP(:,:,i)*DT(:,:,j)' + ...
                      DT(:,:,j)*P*DT(:,:,i)'+DT(:,:,i)*P*DT(:,:,j)'+ ...
                      reshape(D2T(:,jcount),size(T))*P*T'+T*P*reshape(D2T(:,jcount),size(T))' + ...
                      dyn_unvech(D2Om(:,jcount));
                D2P(:,jcount) = dyn_vech(tmp);
                %                 D2P(:,:,i,j) = D2P(:,:,j,i);
            end
        end
    end

end
