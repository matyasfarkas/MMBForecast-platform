function a = vertcat(varargin) % --*-- Unitary tests --*--

%@info:
%! @deftypefn {Function file} {@var{a} =} vertcat (@var{b},@var{c}, ...)
%! @anchor{horzcat}
%! @sp 1
%! Method of the dseries class.
%! @sp 1
%! This method overloads the vertical concatenation operator, so that
%! two (or more) time series objects containing the same variables
%! can be merged using the following syntax:
%!
%!     a = [b; c; d];
%! @sp 2
%! @strong{Inputs}
%! @sp 1
%! @table @ @var
%! @item b
%! Dynare time series object, instantiated by @ref{dseries}.
%! @item c
%! Dynare time series object, instantiated by @ref{dseries}.
%! @end table
%! @sp 2
%! @strong{Outputs}
%! @sp 1
%! @table @var
%! @item a
%! Dynare time series object.
%! @end table
%! @end deftypefn
%@eod:

% Copyright (C) 2013-2017 Dynare Team
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

if nargin==0
    a = DynSeries();
elseif nargin == 1
    a = varargin{1};
elseif nargin>1
    a = varargin{1};
    for i=2:nargin
        a = vertcat_(a,varargin{i});
    end
end

function d = vertcat_(b, c)
d = NaN;
if ~isequal(frequency(b), frequency(c))
    error('dseries::vertcat: Frequencies must be common!')
end
if ~isequal(vobs(b), vobs(c))
    error('dseries::vertcat: Number of variables must be common!')
end
reorder_variables_in_c = false;
if ~isequal(b.name, c.name)
    [t, idx] = ismember(b.name, c.name);
    if all(t)
        reorder_variables_in_c = true;
    else
        error('dseries::vertcat: Variables must be common!')
    end
end
d = b;
if reorder_variables_in_c
    d.data = [b.data; c.data(:,idx)];
else
    d.data = [b.data; c.data];
end
d.dates = [b.dates; c.dates];

%@test:1
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$ B = [transpose(1:10),2*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$ B_name = {'A1';'A2'};
%$
%$ % Define expected results.
%$ e.init = dates(1,1);
%$ e.freq = 1;
%$ e.name = {'A1';'A2'};
%$ e.data = [A;B];
%$
%$ % Instantiate two time series objects.
%$ ts1 = dseries(A,[],A_name,[]);
%$ ts2 = dseries(B,[],B_name,[]);
%$
%$ % Call the tested method.
%$ ts3 = [ts1;ts2];
%$
%$ % Check the results.
%$
%$ t(1) = dassert(ts3.init,e.init);
%$ t(2) = dassert(ts3.freq,e.freq);
%$ t(3) = dassert(ts3.data,e.data);
%$ t(4) = dassert(ts3.name,e.name);
%$ t(5) = dassert(ts3.nobs,20);
%$ T = all(t);
%@eof:1


%@test:2
%$ % Define a data set.
%$ A = [transpose(1:10),2*transpose(1:10)];
%$ B = [transpose(1:10),2*transpose(1:10)];
%$ C = [transpose(1:10),3*transpose(1:10)];
%$
%$ % Define names
%$ A_name = {'A1';'A2'};
%$ B_name = {'A1';'A2'};
%$ C_name = {'A1';'A2'};
%$
%$ % Define expected results.
%$ e.init = dates(1,1);
%$ e.freq = 1;
%$ e.name = {'A1';'A2'};
%$ e.data = [A;B;C];
%$
%$ % Instantiate two time series objects.
%$ ts1 = dseries(A,[],A_name,[]);
%$ ts2 = dseries(B,[],B_name,[]);
%$ ts3 = dseries(C,[],C_name,[]);
%$
%$ % Call the tested method.
%$ ts4 = [ts1; ts2; ts3];
%$
%$ % Check the results.
%$
%$ t(1) = dassert(ts4.init,e.init);
%$ t(2) = dassert(ts4.freq,e.freq);
%$ t(3) = dassert(ts4.data,e.data);
%$ t(4) = dassert(ts4.name,e.name);
%$ t(5) = dassert(ts4.nobs,30);
%$ T = all(t);
%@eof:2

%@test:3
%$ A = dseries([ones(5,1), 2*ones(5,1)],'1938Q4',{'A1', 'A2'});
%$ B = dseries([2*ones(2,1), ones(2,1)],'1945Q3',{'A2', 'A1'});
%$
%$ try
%$    C = [A; B];
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ % Check the results.
%$ if t(1)
%$    t(2) = dassert(C.data,[ones(7,1), 2*ones(7,1)]);
%$ end
%$ T = all(t);
%@eof:3
