function forecast_data=reshape_ascii_forecast_data(endo_nbr, psize, horizon, ascii_data)
% function forecast_data=reshape_ascii_forecast_data(endo_nbr, psize, horizon, ascii_data)
%
% INPUTS
%    endo_nbr:    number of endogenous
%    psize:       number of percentiles
%    horizon:     forecast horizon
%    ascii_data:  data from .out file created by Dan's C code
%
% OUTPUTS
%    forecast_data:    new 3-d array holding data with error bands
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2011-2012 Dynare Team
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

if psize <= 1
    forecast_data = ascii_data;
    return
end

forecast_data = zeros(psize, horizon, endo_nbr);
for i=1:endo_nbr
    for j=1:psize
        forecast_data(j,:,i) = ascii_data(1+horizon*(j-1):horizon*j,i)';
    end
end
end
