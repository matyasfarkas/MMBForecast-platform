function plot_ms_forecast(M_, options_, forecast, figure_name)
% function plot_ms_forecast(M_, options_, forecast, figure_name)
% plots the forecast from the output from a ms-sbvar
%
% INPUTS
%    M_:          (struct)    model structure
%    options_:    (struct)    options
%    forecast:    (matrix)    in the form (percentile x horizon x nvar ), if banded otherwise
%                             ( horizon x nvar )
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

nc = 2;
nr = 2;
nvars = M_.endo_nbr;
endo_names = M_.endo_names;
var_list = endo_names(1:M_.orig_endo_nbr,:);
names = {};
tex_names = {};
m = 1;
for i = 1:size(var_list)
    tmp = strmatch(var_list(i,:),endo_names,'exact');
    if isempty(tmp)
        error([var_list(i,:) ' isn''t and endogenous variable'])
    end
    tex_name = deblank(M_.endo_names_tex(i,:));
    if ~isempty(tex_name)
        names{m} = deblank(var_list(i,:));
        tex_names{m} = tex_name;
        m = m + 1;
    end
end

dims = size(forecast);
if (length(dims) == 2)
    % Point Forecast (horizon x nvars )
    horizon = dims(1);
    num_percentiles = 1;
elseif (length(dims) == 3)
    % Banded Forecast
    horizon = dims(2);
    num_percentiles = dims(1);
else
    error('The impulse response matrix passed to be plotted does not appear to be the correct size');
end

if num_percentiles == 1
    plot_point_forecast(forecast, nvars, nr, nc, var_list, figure_name, ...
                        options_.graph_save_formats, options_.TeX, names, tex_names, ...
                        [options_.ms.output_file_tag filesep 'Output' filesep 'Forecast']);
else
    plot_banded_forecast(forecast, nvars, nr, nc, var_list, num_percentiles, ...
                         figure_name, options_.graph_save_formats, options_.TeX, names, tex_names, ...
                         [options_.ms.output_file_tag filesep 'Output' filesep 'Forecast']);
end

end

function plot_point_forecast(forecast,nvars,nr,nc,endo_names,figure_name,save_graph_formats,TeX,names,tex_names,dirname)
if nvars > nr*nc
    graph_name = ['MS (1) ' figure_name];
    figure('Name', graph_name);
else
    graph_name = figure_name;
    figure('Name', graph_name);
end
m = 1;
n_fig = 1;
for j=1:nvars
    if m > nr*nc
        graph_name = ['MS (' int2str(n_fig) ') ' figure_name];
        dyn_save_graph(dirname,['MS-forecast-' int2str(n_fig)],...
                       save_graph_formats,TeX,names,tex_names,graph_name);
        n_fig =n_fig+1;
        figure('Name', graph_name);
        m = 1;
    end
    subplot(nr,nc,m);
    vn = deblank(endo_names(j,:));
    plot(forecast(:,j))
    title(vn,'Interpreter','none');
    grid on;
    m = m+1;
end
if m > 1
    dyn_save_graph(dirname,['MS-forecast-' int2str(n_fig)],...
                   save_graph_formats,TeX,names,tex_names,graph_name);
end
end

function plot_banded_forecast(forecast,nvars,nr,nc,endo_names,num_percentiles,figure_name,save_graph_formats,TeX,names,tex_names,dirname)
if nvars > nr*nc
    graph_name = ['MS (1) ' figure_name];
    figure('Name', graph_name);
else
    graph_name = figure_name;
    figure('Name', graph_name);
end
m = 1;
n_fig = 1;
for j=1:nvars
    if m > nr*nc
        graph_name = ['MS (' int2str(n_fig) ') ' figure_name];
        dyn_save_graph(dirname,['MS-forecast-' int2str(n_fig)],...
                       save_graph_formats,TeX,names,tex_names,graph_name);
        n_fig =n_fig+1;
        figure('Name',graph_name);
        m = 1;
    end
    subplot(nr,nc,m);
    vn = deblank(endo_names(j,:));
    for k=1:num_percentiles
        if ceil(num_percentiles/2) == k
            plot(forecast(k,:,j),'LineWidth',1.5)
        else
            plot(forecast(k,:,j),'LineWidth',1.1)
        end
        if k==1
            hold on;
        end
    end
    title(vn,'Interpreter','none');
    hold off
    grid on;
    m = m+1;
end
if m > 1
    dyn_save_graph(dirname,['MS-forecast-' int2str(n_fig)],...
                   save_graph_formats,TeX,names,tex_names,graph_name);
end
end
