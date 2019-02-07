%% Function description
% Plot the historical variance decomposition (HVD)
% Shock groups: monetary, demand, supply
% Endogenous variables: GDP, inflation, nominal interest rate

%% Pipeline
% 1. Move to the folder which stores the estimation result
% 2. Identify and group shocks from "M_.exo_names"
% 3. Identify and locate endo variables from "M_.endo_names"
% 4. Retreive the HVD by running "shock_decomposition"
% 5. Save the HVD of the current model to disk
% 6. Create a HVD plot
% 7. Save the HVD plot to disk (if needed)

%% Preambles (same in plot_irf and plot_hvd)

% Find out the model name and vintage date
current_model_info = convertCharsToStrings(M_.fname);
split_position = regexp(current_model_info, "\d{5}");
current_model_name = extractBefore(current_model_info, split_position - 1);
vintage_year = extractBetween(current_model_info, split_position, split_position + 3);
vintage_quarter = extractAfter(current_model_info, split_position + 3);

% Model names
model_names = [...
    "BVAR_MP";
    "BVAR_GLP";
    "US_SW07";
    "US_DNGS14";
    "NK_RW97"];

% Locate current model
current_model_index = find(current_model_name == model_names);

% Model descriptions
model_descriptions = [...
    "BVAR with Minnesta prior";
    "BVAR with Giannone Lenza Primiceri prior";
    "Smets Wouters (2007) US model ";
    "Del Negro Schorfheide Giannoni with credit spread";
    "Small NK US model"];

% Locate current folder (absolute path)
pwd_folder = pwd;

% Current model folder
current_model_folder = pwd_folder + "\..\MODELS\" + current_model_name + "\" + current_model_info;

% Model output folders
model_output_folders = [...
    "BVAR_MP";
    "BVAR_GLP";
    "US_SW07_MHforecast";
    "US_DNGS14_MHforecast";
    "NK_RW97_MHforecast"];
model_output_folders = pwd_folder + "\..\OUTPUT\USMODELS\" + model_output_folders;

% Charts folder
charts_output_folder = pwd_folder + "\..\OUTPUT\Charts";

%% Preambles (different in plot_irf and plot_hvd)

% Length of the HVD
periods_length = 10;
periods = (1:periods_length);

% Headers and subheaders
plot_type = "Historical variance decompositions: ";
headers = plot_type + model_descriptions;
subheaders = ["GDP", "Inflation", "Nominal Interest Rate"];

% Shock names
shock_names = strings(1, size(M_.exo_names, 1));
for shock_index = 1:size(M_.exo_names, 1)
    shock_names(shock_index) = strtrim(convertCharsToStrings(M_.exo_names(shock_index, :)));
end
shock_names = [shock_names, "init"];

% Shock groups
% G1: Supply, G2: Demand, G3: Monetary, G4: Fiscal, G5: Init_condition
% Give the definition of the VD matrix, always end shock_groups with [4, 0]
shock_group_names = ["supply_shocks", "demand_shocks", "monetary_shocks", "initial_conditions"];
switch current_model_name
    case "NK_RW97"
        shock_group_identifiers = [3, 2, 1, 2, 1, 4, 0];
    case "US_SW07"
        shock_group_identifiers = [1, 2, 2, 1, 1, 1, 3, 4, 0];
    case "US_DNGS14"
        shock_group_identifiers = [1, 2, 2, 1, 1, 1, 1, 3, 4, 0];
end

% Variable names and locations
var_names = ["xgdp_q_obs", "pgdp_q_obs", "rff_q_obs"];
var_locations = zeros(1,3);
for all_var_index = 1:size(M_.endo_names, 1)
    all_var_name = convertCharsToStrings(M_.endo_names(all_var_index, :));
    for var_index = 1:3
        if startsWith(all_var_name, var_names(var_index))
            var_locations(var_index) = all_var_index;
        end
    end
end

% Write to filename
write_file_name = "HVD_" + current_model_info + ".xlsx";

%% Run shock_decomposition and retreive HVD
% To run shock_decomposition we need two files in the model folder:
% "current_model_info_dynamic.m" and "current_model_info_stic.m"

% Move to the folder where the results is stored
cd(current_model_folder)

% Retrieve the historical variance decomposition
% hvd_all contains the absolute contributions of ungrouped shocks to three observables
% size(hvd_all) = (n_vars, n_shocks, n_periods)
options_.no_graph.shock_decomposition = 1; % suppress variance decomposition plot from Dynare
hvd_all = shock_decomposition(M_, oo_, options_, '', bayestopt_, estim_params_);
hvd_all = hvd_all.shock_decomposition;
hvd_all = hvd_all(var_locations, :, end:-1:end - periods_length +1);

% hvda contains the absolute contributions of shocks to three observables
% size(hvda) = (n_vars, n_shock_groups, n_periods)
hvda = zeros(length(var_names), length(shock_group_names), periods_length);
for period = periods
    for shock_index = 1:size(hvd_all(:,:,period), 2) - 1
        hvda(:, shock_group_identifiers(shock_index), period) = hvda(:, shock_group_identifiers(shock_index), period) + hvd_all(:, shock_index, period);
    end
end

% size(hvda) = (n_periods, n_shock_groups, n_vars)
hvda_temp = shiftdim(hvda, 1);
hvda = zeros(periods_length, length(shock_group_names), length(var_names));
for var_index = 1:length(var_names)
    hvda(:, :, var_index) = transpose(hvda_temp(:, :, var_index));
end

% hvdr contains the relative contributions of shocks to three observables
% size(hvdr) = size(hvda);
hvdr = zeros(size(hvda));
for var_index = 1:length(var_names)
    hvdr(:, :, var_index) = abs(hvda(:, :, var_index))./sum(abs(hvda(:, :, var_index)), 2);
end

% HVDs contains both relative and absolute variance decompositions
HVDs = struct();
for var_index = 1:length(var_names)
    HVDs.(current_model_name).absolute.(var_names(var_index)) = array2table(hvda(:,:,var_index), "VariableNames", cellstr(shock_group_names));
    HVDs.(current_model_name).relative.(var_names(var_index)) = array2table(hvdr(:,:,var_index), "VariableNames", cellstr(shock_group_names));
end

% Move back to the pwd folder
cd(pwd_folder)

%% Save the HVDs to disk
cd(model_output_folders(current_model_index))
for var_name = var_names
    writetable(HVDs.(current_model_name).absolute.(var_name), write_file_name, "Sheet", "absolute_" + var_name);
    writetable(HVDs.(current_model_name).relative.(var_name), write_file_name, "Sheet", "relative_" + var_name);
end
cd(pwd_folder)

%% Load the HVDs of other models (if exists)
for model_index = 1:5
    if model_index ~= current_model_index
        try
            cd(model_output_folders(model_index))
            read_file_name = "HVD_" + model_names(model_index) + "_" + vintage_year + vintage_quarter + ".xlsx";
            for var_name = var_names
                table_absolute = readtable(read_file_name, "Sheet", "absolute_" + var_name);
                table_relative = readtable(read_file_name, "Sheet", "relative_" + var_name);
                HVDs.(model_names(model_index)).absolute.(var_name) = table_absolute;
                HVDs.(model_names(model_index)).relative.(var_name) = table_relative;
            end
        end
    end
end
cd(pwd_folder)

%% Plot the IRFs of the estimated model
begin_Xs = [-0, -0.175, -0.25, -0.4, -0.4];
diff_Xs = [0, 0.35, 0.25, 0.2, 0.2];

bar_width = 0.6/length(fieldnames(HVDs));
begin_X = begin_Xs(length(fieldnames(HVDs)));
diff_X = diff_Xs(length(fieldnames(HVDs)));

face_alpha = 1;
edge_width = 1.2;
edge_colors = [54, 54, 54; 
    139 0 0;
    0 0 139;
    139 0 139;
    0 139 139]./255;
face_colors = {[255 231 186]/255; [176 224 230]/255; [255 228 225]/255; [207 207 207]/255};

models_for_legend = [];
legend_name_models = strings();

figure;

for current_plot = 1:3
    subplot(2,2,current_plot)
    hold on
    
    % Plot HVDs
    i = 1;
    for model_index = 1:5
        try
            matrix = table2array(HVDs.(model_names(model_index)).absolute.(var_names(current_plot)));
                
            matrix_neg = matrix;
            matrix_neg(matrix_neg > 0) = 0;
            matrix_pos = matrix;
            matrix_pos(matrix_pos < 0) = 0;
            
            hp(i, :) = bar(matrix_pos, bar_width, "stacked");
            hn(i, :) = bar(matrix_neg, bar_width, "stacked");
            
            set(hp(i, :), "XData", 1 + begin_X + diff_X*(i-1):1:periods_length + begin_X + diff_X*(i-1), {"FaceColor"}, face_colors, ...
                "FaceAlpha", face_alpha, "EdgeColor", edge_colors(i, :), "EdgeAlpha", 1, "LineWidth", edge_width);
            set(hn(i, :), "XData", 1 + begin_X + diff_X*(i-1):1:periods_length + begin_X + diff_X*(i-1), {"FaceColor"}, face_colors, ...
                "FaceAlpha", face_alpha, "EdgeColor", edge_colors(i, :), "EdgeAlpha", 1, "LineWidth", edge_width);
                
            models_for_legend(i) = bar(0);
            set(models_for_legend(i), "FaceColor", "w", "EdgeColor", edge_colors(i, :), "EdgeAlpha", 1, "LineWidth", edge_width);
            legend_name_models(i) = model_names(model_index);
                
            i = i + 1;
        end
    end
    
    % Subplot properties
    ax = gca;
    ax.FontSize = 12;
    ax.XLim = [0 periods_length+1];
    % ax.YLim = [0 1];
    ax.XTick = [1:periods_length];
    ax.XGrid = "on";
    ax.YGrid = "on";
    ax.Box = "on";
    ax.Title.String = subheaders(current_plot);
    
end

% Legends for shock groups
subplot(2,2,1)
bars_for_legend = bar(zeros(2, 4), bar_width, "stacked");
set(bars_for_legend, {"FaceColor"}, face_colors, "EdgeColor", "w");   
legend_bars = legend(bars_for_legend, shock_group_names);
set(legend_bars, "Interpreter", "none");
set(legend_bars, "Position", [0.745 0.066 0.16 0.341])
set(legend_bars, "Units","normalized")
title(legend_bars, "Shock groups")

% Legend for models
legend_models = legend(models_for_legend, legend_name_models);
set(legend_models, "Interpreter", "none");
title(legend_models, "Models")

% Adjust position
set(subplot(2,2,1), "Position", [0.131 0.495 0.335 0.341])
set(subplot(2,2,2), "Position", [0.570 0.495 0.335 0.341])
set(subplot(2,2,3), "Position", [0.131 0.066 0.335 0.341])
set(legend_models, "Position", [0.570 0.066 0.16 0.341])
set(gcf, "Units", "normalized", "outerposition", [0 0 1 1])

% Add header
annotation(...
    "textbox", [0.0007 0.891 0.95 0.1], ...
    "String", ["\fontsize{18}", headers(current_model_index), ...
    "\fontsize{12}", "using " + vintage_year + "Q" + vintage_quarter + " vintage"], ...
    "Interpreter", "tex", "EdgeColor", "none", "FontUnits", "normalized", "HorizontalAlignment", "center")

%% Save the IRF plot to disk
choice = questdlg("Would you like to store the graph?", "Save to disk", "Yes", "No", "Yes");
if convertCharsToStrings(choice) == "Yes"
    cd(charts_output_folder)
    saveas(gcf, "HVD_" + current_model_name + "_" + vintage_year + "Q" + vintage_quarter + ".png")
    cd(pwd_folder)
end