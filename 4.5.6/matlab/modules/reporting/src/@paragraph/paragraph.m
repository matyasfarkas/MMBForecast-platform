function o = paragraph(varargin)
%function o = paragraph(varargin)
% Instantiates a paragraph object

% Copyright (C) 2014-2015 Dynare Team
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

o.balancedCols = false;
o.cols = 1;
o.heading = '';
o.indent = true;
o.text = '';

if nargin == 1
    assert(isa(varargin{1}, 'paragraph'),['With one arg to Paragraph constructor, ' ...
                        'you must pass a paragraph object']);
    o = varargin{1};
    return;
elseif nargin > 1
    if round(nargin/2) ~= nargin/2
        error(['@paragraph.paragraph: options must be supplied in name/value ' ...
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
            error('@paragraph.paragraph: %s is not a recognized option.', pair{1});
        end
    end
end

assert(islogical(o.indent), '@paragraph.paragraph: indent must be either true or false');
assert(islogical(o.balancedCols), '@paragraph.paragraph: balancedCols must be either true or false');
assert(isint(o.cols), '@paragraph.paragraph: cols must be an integer');
assert(ischar(o.text), '@paragraph.paragraph: text must be a string');
assert(ischar(o.heading), '@paragraph.paragraph: heading must be a string');

% Create paragraph object
o = class(o, 'paragraph');
end
