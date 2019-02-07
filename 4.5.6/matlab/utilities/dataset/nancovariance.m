function CovarianceMatrix = nancovariance(data)
% Computes the covariance matrix of a sample (possibly with missing observations).

%@info:
%! @deftypefn {Function File} {@var{CovarianceMatrix} =} compute_corr(@var{data})
%! @anchor{compute_cova}
%! This function computes covariance matrix of a sample defined by a dseries object (possibly with missing observations).
%!
%! @strong{Inputs}
%! @table @var
%! @item data
%! a T*N array of real numbers.
%! @end table
%!
%! @strong{Outputs}
%! @table @var
%! @item CovarianceMatrix
%! Array of real numbers.
%! @end table
%!
%! @strong{This function is called by:}
%! @ref{descriptive_statistics}.
%!
%! @strong{This function calls:}
%! @ref{ndim}, @ref{demean}, @ref{nandemean}.
%!
%! @strong{Remark 1.} On exit, a new field is appended to the structure: @code{dataset_.descriptive.cova} is a
%! @tex{n\times n} vector (where @tex{n} is the number of observed variables as defined by @code{dataset_.info.nvobs}).
%!
%! @end deftypefn
%@eod:

% Copyright (C) 2011-2014 Dynare Team
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

% Initialize the output.
CovarianceMatrix = zeros(size(data,2));

if isanynan(data)
    data = bsxfun(@minus,data,nanmean(data));
    for i=1:size(data,2)
        for j=i:size(data,2)
            CovarianceMatrix(i,j) = nanmean(data(:,i).*data(:,j));
            if j>i
                CovarianceMatrix(j,i) = CovarianceMatrix(i,j);
            end
        end
    end
else
    data = bsxfun(@minus,data,mean(data));
    CovarianceMatrix = (transpose(data)*data)/size(data,1);
end

%@test:1
%$
%$ % Define a dataset.
%$ data1 = randn(10000000,2);
%$
%$ % Same dataset with missing observations.
%$ data2 = data1;
%$ data2(45,1) = NaN;
%$ data2(57,2) = NaN;
%$ data2(367,:) = NaN(1,2);
%$
%$ t = zeros(2,1);
%$
%$ % Call the tested routine.
%$ try
%$     c1 = nancovariance(data1);
%$     t(1) = 1;
%$ catch
%$     t(1) = 0;
%$ end
%$ try
%$     c2 = nancovariance(data2);
%$     t(2) = 1;
%$ catch
%$     t(2) = 0;
%$ end
%$
%$ if t(1) && t(2)
%$     t(3) = max(max(abs(c1-c2)))<1e-4;
%$ end
%$
%$ % Check the results.
%$ T = all(t);
%@eof:1