function o = page(varargin)
%function o = page(varargin)
% Page Class Constructor
%
% INPUTS
%   0 args => empty page
%   1 arg (page class) => copy object
%
% OUTPUTS
%   none
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
o.paper = '';
o.title = {''};
titleFormatDefalut = {'\large\bfseries'};
o.titleFormat = titleFormatDefalut;
o.titleTruncate = '';
o.orientation = '';
o.footnote = {};
o.sections = {};

o.pageDirName = 'tmpRepDir';
o.latex = '';

if nargin == 1
    assert(isa(varargin{1}, 'page'), ['@page.page: with one arg to Page ' ...
                        'constructor, you must pass a page object']);
    o = varargin{1};
    return;
elseif nargin > 1
    if round(nargin/2) ~= nargin/2
        error(['@page.page: options must be supplied in name/value ' ...
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
            error('@page.page: %s is not a recognized option.', pair{1});
        end
    end
end

% Check options provided by user
if ischar(o.title)
    o.title = {o.title};
end
if ischar(o.titleFormat)
    o.titleFormat = {o.titleFormat};
end
if length(o.title) ~= length(o.titleFormat)
    o.titleFormat = repmat(titleFormatDefalut, 1, length(o.title));
end
assert(iscellstr(o.title), ...
       '@page.page: title must be a cell array of strings');
assert(iscellstr(o.titleFormat), ...
       '@page.page: titleFormat must be a cell array of strings');
assert((ischar(o.titleTruncate) && isempty(o.titleTruncate)) || ...
       isint(o.titleTruncate), ...
       '@page.page: titleTruncate must be empty or an integer.');
assert(ischar(o.pageDirName), '@page.page: pageDirName must be a string');
assert(ischar(o.latex), ...
       '@page.page: latex must be a string');
valid_paper = {'a4', 'letter'};
assert(any(strcmp(o.paper, valid_paper)), ...
       ['@page.page: paper must be one of ' addCommasToCellStr(valid_paper)]);

valid_orientation = {'portrait', 'landscape'};
assert(any(strcmp(o.orientation, valid_orientation)), ...
       ['@page.page: orientation must be one of ' addCommasToCellStr(valid_orientation)]);

if ischar(o.footnote)
    o.footnote = {o.footnote};
end
assert(iscellstr(o.footnote), ...
       '@page.page: footnote must be a cell array of string(s)');

% Create page object
o = class(o, 'page');
end
