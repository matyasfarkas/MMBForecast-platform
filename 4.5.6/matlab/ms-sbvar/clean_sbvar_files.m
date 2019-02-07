function clean_sbvar_files()
% function clean_sbvar_files()
% Remove files created by sbvar
%
% INPUTS
%   none
%
% OUTPUTS
%   none
%
% SPECIAL REQUIREMENTS
%   none

% Copyright (C) 2010-2011 Dynare Team
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

delete_if_exists('outm.mat');
delete_if_exists('g1.mat');
delete_if_exists('H.dat');
delete_if_exists('outactqmygdata.prn');
delete_if_exists('outdata_a0dp_const.mat');
delete_if_exists('outyrqm.prn');
end
