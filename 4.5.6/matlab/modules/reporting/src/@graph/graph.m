function o = graph(varargin)
%function o = graph(varargin)
% Graph Class Constructor
%
% INPUTS
%   varargin        0 args  : empty graph object
%                   1 arg   : must be graph object (return a copy of arg)
%                   > 1 args: option/value pairs (see structure below for
%                   options)
%
% OUTPUTS
%   o   [graph] graph object
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

o.series = {};

o.title = '';
o.titleFormat = '';
o.titleFontSize = 'normalsize';
o.ylabel = '';
o.xlabel = '';

o.axisShape = 'box';

o.graphDirName = 'tmpRepDir';
o.graphName = '';
o.data = '';
o.seriesToUse = '';
o.xrange = '';
o.xAxisTight = true;
o.yrange = '';
o.yAxisTight = false;

o.shade = '';
o.shadeColor = 'green';
o.shadeOpacity = 20;

o.showGrid = true;

o.showLegend = false;
o.legendAt = [];
o.showLegendBox = false;
o.legendLocation = 'south east';
o.legendOrientation = 'horizontal';
o.legendFontSize = 'tiny';

o.showZeroline = false;
o.zeroLineColor = 'black';

o.xTicks = [];
o.xTickLabels = {};
o.xTickLabelRotation = 0;
o.xTickLabelAnchor = 'east';

o.yTickLabelScaled = true;
o.yTickLabelPrecision = 0;
o.yTickLabelFixed = true;
o.yTickLabelZeroFill = true;

o.tickFontSize = 'normalsize';

o.width = 6;
o.height = 4.5;

o.miscTikzPictureOptions = '';
o.miscTikzAxisOptions = '';

o.writeCSV = false;

if nargin == 1
    assert(isa(varargin{1}, 'graph'),['@graph.graph: with one arg you ' ...
                        'must pass a graph object']);
    o = varargin{1};
    return;
elseif nargin > 1
    if round(nargin/2) ~= nargin/2
        error(['@graph.graph: options must be supplied in name/value ' ...
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
            error('@graph.graph: %s is not a recognized option.', pair{1});
        end
    end
end

% Check options provided by user
if ischar(o.title)
    o.title = {o.title};
end
assert(iscellstr(o.title), '@graph.graph: title must be a cell array of string(s)');
assert(ischar(o.titleFormat), '@graph.graph: titleFormat file must be a string');
assert(ischar(o.xlabel), '@graph.graph: xlabel file must be a string');
assert(ischar(o.ylabel), '@graph.graph: ylabel file must be a string');
assert(ischar(o.miscTikzPictureOptions), '@graph.graph: miscTikzPictureOptions file must be a string');
assert(ischar(o.miscTikzAxisOptions), '@graph.graph: miscTikzAxisOptions file must be a string');
assert(ischar(o.graphName), '@graph.graph: graphName must be a string');
assert(ischar(o.graphDirName), '@graph.graph: graphDirName must be a string');
assert(islogical(o.showGrid), '@graph.graph: showGrid must be either true or false');
assert(islogical(o.xAxisTight), '@graph.graph: xAxisTight must be either true or false');
assert(islogical(o.yAxisTight), '@graph.graph: yAxisTight must be either true or false');
assert(islogical(o.showLegend), '@graph.graph: showLegend must be either true or false');
assert(isempty(o.legendAt) || (isfloat(o.legendAt) && length(o.legendAt)==2), ...
       '@graph.graph: legendAt must be a double array of size two');
assert(islogical(o.showLegendBox), '@graph.graph: showLegendBox must be either true or false');
assert(islogical(o.showZeroline), '@graph.graph: showZeroline must be either true or false');
assert(isfloat(o.shadeOpacity) && length(o.shadeOpacity)==1 && ...
       o.shadeOpacity >= 0 && o.shadeOpacity <= 100, ...
       '@graph.graph: o.shadeOpacity must be a real in [0 100]');
assert(isfloat(o.width), '@graph.graph: o.width must be a real number');
assert(isfloat(o.height), '@graph.graph: o.height must be a real number');
assert(isfloat(o.xTickLabelRotation), '@graph.graph: o.xTickLabelRotation must be a real number');
assert(ischar(o.xTickLabelAnchor), '@graph.graph: xTickLabelAnchor must be a string');
assert(isint(o.yTickLabelPrecision), '@graph.graph: o.yTickLabelPrecision must be an integer');
assert(islogical(o.yTickLabelFixed), '@graph.graph: yTickLabelFixed must be either true or false');
assert(islogical(o.yTickLabelZeroFill), '@graph.graph: yTickLabelZeroFill must be either true or false');
assert(islogical(o.yTickLabelScaled), '@graph.graph: yTickLabelScaled must be either true or false');
assert(islogical(o.writeCSV), '@graph.graph: writeCSV must be either true or false');
assert(ischar(o.shadeColor), '@graph.graph: shadeColor must be a string');
assert(ischar(o.zeroLineColor), '@graph.graph: zeroLineColor must be a string');
assert(any(strcmp(o.axisShape, {'box', 'L'})), ['@graph.graph: axisShape ' ...
                    'must be one of ''box'' or ''L''']);
valid_legend_locations = ...
    {'south west','south east','north west','north east','outer north east'};
assert(any(strcmp(o.legendLocation, valid_legend_locations)), ...
       ['@graph.graph: legendLocation must be one of ' addCommasToCellStr(valid_legend_locations)]);

valid_font_sizes = {'tiny', 'scriptsize', 'footnotesize', 'small', ...
                    'normalsize', 'large', 'Large', 'LARGE', 'huge', 'Huge'};
assert(any(strcmp(o.legendFontSize, valid_font_sizes)), ...
       ['@graph.graph: legendFontSize must be one of ' addCommasToCellStr(valid_font_sizes)]);
assert(any(strcmp(o.titleFontSize, valid_font_sizes)), ...
       ['@graph.graph: titleFontSize must be one of ' addCommasToCellStr(valid_font_sizes)]);
assert(any(strcmp(o.tickFontSize, valid_font_sizes)), ...
       ['@graph.graph: tickFontSize must be one of ' addCommasToCellStr(valid_font_sizes)]);

valid_legend_orientations = {'vertical', 'horizontal'};
assert(any(strcmp(o.legendOrientation, valid_legend_orientations)), ...
       ['@graph.graph: legendOrientation must be one of ' addCommasToCellStr(valid_legend_orientations)]);

assert(isempty(o.shade) || (isdates(o.shade) && o.shade.ndat >= 2), ...
       ['@graph.graph: shade is specified as a dates range, e.g. ' ...
        '''dates(''1999q1''):dates(''1999q3'')''.']);
assert(isempty(o.xrange) || (isdates(o.xrange) && o.xrange.ndat >= 2), ...
       ['@graph.graph: xrange is specified as a dates range, e.g. ' ...
        '''dates(''1999q1''):dates(''1999q3'')''.']);
assert(isempty(o.yrange) || (isfloat(o.yrange) && length(o.yrange) == 2 && ...
                             o.yrange(1) < o.yrange(2)), ...
       ['@graph.graph: yrange is specified an array with two float entries, ' ...
        'the lower bound and upper bound.']);
assert(isempty(o.data) || isdseries(o.data), ['@graph.graph: data must ' ...
                    'be a dseries']);
assert(isempty(o.seriesToUse) || iscellstr(o.seriesToUse), ['@graph.graph: ' ...
                    'seriesToUse must be a cell array of string(s)']);
assert(isempty(o.xTicks) || isfloat(o.xTicks),...
       '@graph.graph: xTicks must be a numerical array');
assert(iscellstr(o.xTickLabels) || (ischar(o.xTickLabels) && strcmpi(o.xTickLabels, 'ALL')), ...
       ['@graph.graph: xTickLabels must be a cell array of strings or ' ...
        'equivalent to the string ''ALL''']);
if ~isempty(o.xTickLabels)
    assert((ischar(o.xTickLabels) && strcmpi(o.xTickLabels, 'ALL')) || ...
           ~isempty(o.xTicks), ['@graph.graph: if you set xTickLabels and ' ...
                        'it''s not equal to ''ALL'', you must set xTicks']);
end
if ~isempty(o.xTicks)
    assert(~isempty(o.xTickLabels), '@graph.graph: if you set xTicks, you must set xTickLabels');
end

% using o.seriesToUse, create series objects and put them in o.series
if ~isempty(o.data)
    if isempty(o.seriesToUse)
        for i=1:o.data.vobs
            o.series{end+1} = report_series('data', o.data{o.data.name{i}});
        end
    else
        for i=1:length(o.seriesToUse)
            o.series{end+1} = report_series('data', o.data{o.seriesToUse{i}});
        end
    end
end
o = rmfield(o, 'seriesToUse');
o = rmfield(o, 'data');

if ~exist(o.graphDirName, 'file')
    mkdir(o.graphDirName);
end

% Create graph object
o = class(o, 'graph');
end
