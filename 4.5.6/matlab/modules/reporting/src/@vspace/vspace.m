function o = vspace(varargin)
%function o = vspace(varargin)
% Vspace Class Constructor
%
% INPUTS
%   0 args => empty vspace
%   1 arg (vspace class) => copy object
%
% OUTPUTS
%   none
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

o = struct;
o.number = 1;
o.hline = 0;

if nargin == 1
    assert(isa(varargin{1}, 'vspace'), ['@vspace.vspace: with one arg to Vspace ' ...
                        'constructor, you must pass a vspace object']);
    o = varargin{1};
    return;
elseif nargin > 1
    if round(nargin/2) ~= nargin/2
        error(['@vspace.vspace: options must be supplied in name/value ' ...
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
            error('@vspace.vspace: %s is not a recognized option.', pair{1});
        end
    end
end

% Check options provided by user
assert(isint(o.number), '@vspace.vspace: number must be an integer');
assert(isint(o.hline), '@vspace.vspace: hline must be an integer');

% Create vspace object
o = class(o, 'vspace');
end
