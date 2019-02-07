function plot_ms_probabilities(computed_probabilities, options_)
% function plot_ms_probabilities(computed_probabilities, options_)
% Plots the regime probablities for each graph
%
% INPUTS
%    computed_probabilities:      Txnstates
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

[T,num_grand_regimes] = size(computed_probabilities);
num_chains = length(options_.ms.ms_chain);
for i=1:num_chains
    chains(i).num_regimes = length(options_.ms.ms_chain(i).regime);
    chains(i).probabilities = zeros([T,chains(i).num_regimes]);
end

for t=1:T
    chains = iterate_chain(computed_probabilities(t,:), t, chains, 1, num_chains);
end

for i=1:num_chains
    graph_name = ['MS-Probabilities, Chain ' int2str(i)];
    figure('Name',graph_name)
    plot(chains(i).probabilities,'LineWidth', 1.2);
    ltxt = {};
    for j=1:chains(i).num_regimes
        ltxt{j} = ['Regime ' int2str(j)];
    end
    legend(ltxt{:});
    title(['Chain ' int2str(i)]);
    ylim([0 1.0]);
    dyn_save_graph([options_.ms.output_file_tag filesep 'Output' filesep 'Probabilities'], ['MS-Probabilities-Chain-' int2str(i)], ...
                   options_.graph_save_formats,options_.TeX,[],[],graph_name);
end
end

function [chains] = iterate_chain(probs, t, chains, chain, num_chains)
offset_length = length(probs)/chains(chain).num_regimes;
for i=1:chains(chain).num_regimes
    p = probs( (i-1)*offset_length+1 : i*offset_length );
    chains(chain).probabilities(t, i) = chains(chain).probabilities(t, i) + sum( p );
    if chain < num_chains
        chains = iterate_chain(p, t, chains, chain+1, num_chains);
    end
end
end
