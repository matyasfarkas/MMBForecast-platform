function covM = fn_construct_cov(diagcovM, corM)
%  function corM = corr(covM)
%  diagcovM: diag of covariance matrix (input)
%  corM:  correlation matrix  (input)
%  covM:  covariance matrix (output)
%  9 May 2008
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


[nvar,jnk]=size(corM);
if nvar~=jnk
	error('The input matrix must be square!!')
end
covM=zeros(nvar);

for i=1:nvar
   for j=1:nvar
     covM(i,j)= corM(i,j) * sqrt( diagcovM(i) * diagcovM(j) );
   end
end

