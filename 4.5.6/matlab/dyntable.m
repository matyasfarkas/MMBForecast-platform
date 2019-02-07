function dyntable(options_,title,headers,labels,values,label_width,val_width,val_precis,optional_header)
% function dyntable(title,headers,labels,values,label_width,val_width,val_precis)
% Inputs:
%   options_    [structure]         Dynare options structure
%   title       [string]            Table title
%   headers     [n by nchar]        character array of labels for header row
%   labels      [n by nchar]        character array of labels for label column
%   values      [matrix]            matrix of values to display
%   label_width [scalar]            Width of the label
%   val_width   [scalar]            Width of value column
%   val_precis  [integer]           precision of displayed values
%
%
% Copyright (C) 2002-2017 Dynare Team
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

if options_.noprint
    return
end

%% get width of label column
if ~isempty(label_width)
    label_width = max(size(deblank(char(headers(1,:),labels)),2)+2, ...
                      label_width);
else %use default length
    label_width = max(size(deblank(char(headers(1,:),labels)),2))+2;
end
label_format_leftbound  = sprintf('%%-%ds',label_width);

%% get width of label column
if all(~isfinite(values))
    values_length = 4;
else
    values_length = max(ceil(max(max(log10(abs(values(isfinite(values))))))),1)+val_precis+1;
end
if any(values) < 0 %add one character for minus sign
    values_length = values_length+1;
end

%% get width of header strings
headers_length = max(size(deblank(headers(2:end,:)),2));
if ~isempty(val_width)
    val_width = max(max(headers_length,values_length)+2,val_width);
else
    val_width = max(headers_length,values_length)+2;
end
value_format  = sprintf('%%%d.%df',val_width,val_precis);
header_string_format  = sprintf('%%%ds',val_width);

if length(title) > 0
    fprintf('\n\n%s\n',title);
end
%Create and print header string
if nargin==9
    disp(optional_header)
end
if length(headers) > 0
    hh = sprintf(label_format_leftbound ,deblank(headers(1,:)));
    for i=2:size(headers,1)
        hh  = [hh sprintf(header_string_format,deblank(headers(i,:)))];
    end
    disp(hh);
end

for i=1:size(values,1)
    disp([sprintf(label_format_leftbound ,deblank(labels(i,:))) sprintf(value_format ,values(i,:))]);
end