function o = simpsaget(options,name,default,flag)
%SIMPSAGET Get SIMPSA OPTIONS parameters.
%   VAL = SIMPSAGET(OPTIONS,'NAME') extracts the value of the named parameter
%   from optimization options structure OPTIONS, returning an empty matrix if
%   the parameter value is not specified in OPTIONS.  It is sufficient to
%   type only the leading characters that uniquely identify the
%   parameter.  Case is ignored for parameter names.  [] is a valid OPTIONS
%   argument.
%
%   VAL = SIMPSAGET(OPTIONS,'NAME',DEFAULT) extracts the named parameter as
%   above, but returns DEFAULT if the named parameter is not specified (is [])
%   in OPTIONS.  For example
%
%     val = simpsaget(opts,'TolX',1e-4);
%
%   returns val = 1e-4 if the TolX property is not specified in opts.
%
%   See also SIMPSASET, SIMPSA

% Copyright (C) 2006 Brecht Donckels, BIOMATH, brecht.donckels@ugent.be
% Copyright (C) 2013-2017 Dynare Team.
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


% undocumented usage for fast access with no error checking
if (nargin == 4) && isequal(flag,'fast')
    o = getknownfield(options,name,default);
    return
end

if nargin < 2
    error('MATLAB:odeget:NotEnoughInputs','Not enough input arguments.');
end
if nargin < 3
    default = [];
end

if ~isempty(options) && ~isa(options,'struct')
    error('MATLAB:odeget:Arg1NotODESETstruct',...
          'First argument must be an options structure created with ODESET.');
end

if isempty(options)
    o = default;
    return
end

Names = [
    'TEMP_START               '
    'TEMP_END                 '
    'COOL_RATE                '
    'INITIAL_ACCEPTANCE_RATIO '
    'MIN_COOLING_FACTOR       '
    'MAX_ITER_TEMP_FIRST      '
    'MAX_ITER_TEMP_LAST       '
    'MAX_ITER_TEMP            '
    'MAX_ITER_TOTAL           '
    'MAX_TIME                 '
    'MAX_FUN_EVALS            '
    'TOLX                     '
    'TOLFUN                   '
    'DISPLAY                  '
    'OUTPUT_FCN               '
        ];

names = lower(Names);

lowName = lower(name);
j = strmatch(lowName,names);
if isempty(j)               % if no matches
    error('MATLAB:odeget:InvalidPropName',...
          ['Unrecognized property name ''%s''.  ' ...
           'See ODESET for possibilities.'], name);
elseif length(j) > 1            % if more than one match
                                % Check for any exact matches (in case any names are subsets of others)
    k = strmatch(lowName,names,'exact');
    if length(k) == 1
        j = k;
    else
        msg = sprintf('Ambiguous property name ''%s'' ', name);
        msg = [msg '(' deblank(Names(j(1),:))];
        for k = j(2:length(j))'
            msg = [msg ', ' deblank(Names(k,:))];
        end
        msg = sprintf('%s).', msg);
        error('MATLAB:odeget:AmbiguousPropName', msg);
    end
end

if any(strcmp(fieldnames(options),deblank(Names(j,:))))
    o = options.(deblank(Names(j,:)));
    if isempty(o)
        o = default;
    end
else
    o = default;
end

% --------------------------------------------------------------------------
function v = getknownfield(s, f, d)
%GETKNOWNFIELD  Get field f from struct s, or else yield default d.

if isfield(s,f)   % s could be empty.
    v = subsref(s, struct('type','.','subs',f));
    if isempty(v)
        v = d;
    end
else
    v = d;
end
