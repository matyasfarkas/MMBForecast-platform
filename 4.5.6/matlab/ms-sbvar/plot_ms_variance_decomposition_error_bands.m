function plot_ms_variance_decomposition_error_bands(M_, options_, vddata, figure_name)
% function plot_ms_variance_decomposition_error_bands(M_, options_, vddata, figure_name)
% plots the variance decomposition with percentiles
%
% INPUTS
%    M_:          (struct)    model structure
%    options_:    (struct)    options
%    vddata:      (matrix)    variance_decomposition (percentile, options_.ms.horizon, nvar
%                             x nvar)
%    figure_name: (string)    title
%
% OUTPUTS
%    none
%
% SPECIAL REQUIREMENTS
%    none

% Copyright (C) 2011-2017 Dynare Team
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

nvars = M_.endo_nbr;
endo_names = M_.endo_names;
var_list = endo_names(1:M_.orig_endo_nbr,:);

names = {};
tex_names = {};
m = 1;
for i = 1:size(var_list)
    tmp = strmatch(var_list(i,:), endo_names, 'exact');
    if isempty(tmp)
        error([var_list(i,:) ' isn''t an endogenous variable'])
    end
    tex_name = deblank(M_.endo_names_tex(tmp,:));
    if ~isempty(tex_name)
        names{m} = deblank(var_list(i,:));
        tex_names{m} = tex_name;
        m = m + 1;
    end
end

for i=1:M_.exo_nbr
    tex_name = deblank(M_.exo_names_tex(i,:));
    if ~isempty(tex_name)
        names{m} = deblank(M_.exo_names(i,:));
        tex_names{m} = tex_name;
        m = m + 1;
    end
end

dims = size(vddata);
if length(dims) ~= 3
    error('The variance decomposition matrix passed to be plotted does not appear to be the correct size');
end
num_percentiles = dims(1);

if size(endo_names, 1) ~= nvars
    error('The names passed are not the same length as the number of variables')
end

for s=1:nvars
    shock = zeros(options_.ms.horizon, nvars, num_percentiles);
    for n=1:num_percentiles
        for i=1:nvars
            shock(:,i,n) = vddata(n, :, ((i-1) + ((s-1)*nvars)+1));
        end
    end
    plot_banded_vddata_for_shock(shock, nvars, endo_names, ...
                                 deblank(endo_names(s,:)), figure_name, ...
                                 [options_.ms.output_file_tag filesep 'Output' filesep 'Variance_Decomposition'], ...
                                 options_, names, tex_names);
end
end

function [fig] = plot_banded_vddata_for_shock(vddata, nvars, endo_names, ...
                                              shock_name, figure_name, dirname, options_, names, tex_names)
fig = figure('Name', figure_name);
npercentiles = size(vddata,3);
for k=1:nvars
    subplot(ceil(sqrt(nvars)), ceil(sqrt(nvars)),k);
    for nn=1:npercentiles
        plot(vddata(:,k,nn))
        hold on
    end
    hold off
    disp([endo_names(k,:) ' contribution to ' shock_name]);
    title([endo_names(k,:) ' contribution to ' shock_name]);
end
dyn_save_graph(dirname, [figure_name ' ' shock_name], ...
               options_.graph_save_formats, options_.TeX, names, tex_names, ...
               [figure_name ' ' shock_name]);
end
