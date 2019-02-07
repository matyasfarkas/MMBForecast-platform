function h = dyn_figure(nodisplay, varargin)
%function h = dyn_figure(nodisplay, varargin)
% initializes figures for DYNARE
%
% INPUTS
%    nodisplay: the value of the command-specific nodisplay argument or options_.nodisplay
%    varargin: the same list of possible inputs of the MATLAB function figure
%
% OUTPUTS
%    h     : figure handle
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2012-2017 Dynare Team
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

if nodisplay
    h = figure(varargin{:},'visible','off');
else
    h = figure(varargin{:});
end
