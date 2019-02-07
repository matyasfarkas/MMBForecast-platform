function o = report_series(varargin)
%function o = report_series(varargin)
% Report_Series Class Constructor
%
% INPUTS
%   varargin        0 args  : empty report_series object
%                   1 arg   : must be report_series object (return a copy of arg)
%                   > 1 args: option/value pairs (see structure below for
%                   options)
%
% OUTPUTS
%   o   [report_series] report_series object
%
% SPECIAL REQUIREMENTS
%   none

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

o = struct;

o.data = '';

o.graphFanShadeColor = '';
o.graphFanShadeOpacity = 50;

o.graphLegendName = '';

o.graphLineColor = 'black';
o.graphLineStyle = 'solid';
o.graphLineWidth = 0.5;

o.graphShowInLegend = true;

o.graphMarker = '';
o.graphMarkerEdgeColor = '';
o.graphMarkerFaceColor = '';
o.graphMarkerSize = 1;

o.graphMiscTikzAddPlotOptions = '';

o.graphHline = {};
o.graphVline = dates();

o.graphBar = false;
o.graphBarColor = 'black';
o.graphBarFillColor = 'black';
o.graphBarWidth = 2;

o.tableShowMarkers = false;
o.tableNegColor = 'red';
o.tablePosColor = 'blue';
o.tableMarkerLimit = 1e-4;

o.tableSubSectionHeader = '';
o.tableAlignRight = false;

o.tableRowColor = 'white';
o.tableRowIndent = 0;

o.tableDataRhs = '';
o.tableNaNSymb = 'NaN';

o.tablePrecision = '';

o.zeroTol = 1e-6;

if nargin == 1
    assert(isa(varargin{1}, 'report_series'),['@report_series.report_series: with one arg you ' ...
                        'must pass a report_series object']);
    o = varargin{1};
    return;
elseif nargin > 1
    if round(nargin/2) ~= nargin/2
        error(['@report_series.report_series: options must be supplied in name/value ' ...
               'pairs.']);
    end

    optNames = fieldnames(o);

    % overwrite default values
    for pair = reshape(varargin, 2, [])
        ind = find(strcmpi(optNames, pair{1}));
        assert(isempty(ind) || length(ind) == 1);
        if ~isempty(ind)
            o.(optNames{ind}) = pair{2};
        else
            error('@report_series.report_series: %s is not a recognized option.', pair{1});
        end
    end
end

if ~isempty(o.graphLegendName)
    o.data = o.data.tex_rename(o.graphLegendName);
end

% Create report_series object
o = class(o, 'report_series');
end
