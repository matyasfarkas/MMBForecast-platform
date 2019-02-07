function e = ep_accuracy_check(M,options,oo)

endo_simul = oo.endo_simul;
n = size(endo_simul,2);
[initialconditions, innovations, pfm, ep, verbosity, options, oo] = ...
    extended_path_initialization([], n-1, [], options, M, oo);

options.ep.accuracy.stochastic.order = options.ep.stochastic.order;
[nodes,weights,nnodes] = setup_integration_nodes(options.ep.accuracy,pfm);

e = zeros(M.endo_nbr,n);
for i=1:n
    e(:,i) = euler_equation_error(endo_simul(:,i),oo.exo_simul, ...
                                  innovations, M, options,oo,pfm,nodes,weights);
end