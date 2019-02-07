function o = addSeries(o, varargin)
%function o = addSeries(o, varargin)
% Add a graph to the current section of the current page in the report
%
% INPUTS
%   o          [report]  report object
%   varargin             arguments to @graph/addSeries
%
% OUTPUTS
%   o          [report]  updated report object
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2013-2015 Dynare Team
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

assert(length(o.pages) > 0, ...
       ['@report.addSeries: Before adding a series, you must add a page, ' ...
        'section, and either a graph or a table.']);
assert(length(o.pages{end}.sections) > 0, ...
       ['@report.addSeries: Before adding a series, you must add a section and ' ...
        'either a graph or a table']);
assert(length(o.pages{end}.sections{end}.elements) > 0, ...
       ['@report.addSeries: Before adding a series, you must add ' ...
        'either a graph or a table']);
assert(isa(o.pages{end}.sections{end}.elements{end}, 'graph') || ...
       isa(o.pages{end}.sections{end}.elements{end}, 'report_table'), ...
       '@report.addSeries: you can only add a series to a report_table or graph object');

o.pages{end}.sections{end}.elements{end} = ...
    o.pages{end}.sections{end}.elements{end}.addSeries(varargin{:});
end
