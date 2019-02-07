function title=add_filter_subtitle(title,options_)

% Copyright (C) 2015-2017 Dynare Team
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

if ~options_.hp_filter && ~options_.one_sided_hp_filter  && ~options_.bandpass.indicator %do not filter
                                                                                         %nothing to add here
elseif ~options_.hp_filter && ~options_.one_sided_hp_filter && options_.bandpass.indicator
    title = [title ' (Bandpass filter, (' ...
             num2str(options_.bandpass.passband(1)),' ',num2str(options_.bandpass.passband(2)), '))'];
elseif options_.hp_filter && ~options_.one_sided_hp_filter  && ~options_.bandpass.indicator %filter with HP-filter
    title = [title ' (HP filter, lambda = ' ...
             num2str(options_.hp_filter) ')'];
elseif ~options_.hp_filter && options_.one_sided_hp_filter && ~options_.bandpass.indicator
    title = [title ' (One-sided HP filter, lambda = ' ...
             num2str(options_.one_sided_hp_filter) ')'];
end
end