function rec_dates = fn_ReadRecessionDates57_01()
% Recession dates.  For some reasons, Matlab code will only take pairs of rec_dates when using rectangle.
% Note that -1 in the following is for the graph lie in the correct position.
%   For example, 1961+(1-1)/12 means Jan 1961 plotted right at 1961.
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

rec_dates=[
    1957+(7-1)/12 1958+(3-1)/12
    1960+(3-1)/12 1961+(1-1)/12
    1969+(11-1)/12 1970+(10-1)/12
    1973+(10-1)/12 1975+(2-1)/12
    1980+(1-1)/12 1980+(6-1)/12
    1981+(6-1)/12 1982+(10-1)/12
    1990+(6-1)/12 1991+(2-1)/12
    2001+(2-1)/12 2001+(10-1)/12
    2007+(12-1)/12 2008+(12-1)/12];

