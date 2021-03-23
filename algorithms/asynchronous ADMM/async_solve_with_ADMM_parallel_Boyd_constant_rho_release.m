function [cpx,central_soln] = async_solve_with_ADMM_parallel_Boyd_constant_rho_release(select, A, lhs, rhs, ub, lb, ...
    f, cluster_id_from, cluster_id_to)

cluster_id_from = cluster_id_from';
cluster_id_to = cluster_id_to';
cluster_annotation = full([cluster_id_from'; cluster_id_to']);
row_size = size(A,1);
column_size = size(A,2);

%Zentrales Problem aufbauen 
orig_model = build_central(select, A, lhs, rhs, ub, lb, f);

%(optional) Zentrales Problem l�sen zum Vergleich der Kostenfunktionen
[central_soln, central_soln_obj] = linprog(orig_model.Model.obj, ...
                                           orig_model.Model.A, ...
                                           orig_model.Model.rhs, ...
                                           [], ...
                                           [], ...
                                           orig_model.Model.lb, ...
                                           orig_model.Model.ub);
                                       
%ADMM-Parameter von select lesen
maxiter = select.opt.decomposition.maxiter; 
mu = select.opt.decomposition.mu;
tau = select.opt.decomposition.tau;
Plink_init = select.opt.decomposition.Plink_init;
lambdas_init = select.opt.decomposition.lambdas_init;
rhos_init = select.opt.decomposition.rhos_init;
convergence_tolerance = select.opt.decomposition.convergence_tolerance;
convergence_type = select.opt.decomposition.convergence_type;

tau = 0;

%Ein Array von eindeutigen Clusternamen
cluster_names = unique(cluster_id_from);

%Verschiedene Indizmengen initializieren bzw. definieren
%f�r jeweilige Cluster
global_var_set = [];
global_int_set = [];
internal_var_indices = cell(numel(cluster_names),1);
internal_plus_neighbor_var_indices = cell(numel(cluster_names),1);
coupling_var_indices = cell(numel(cluster_names),1);
relative_indices = cell(numel(cluster_names),1);
local_position_coupling_indices = cell(numel(cluster_names),1);

for cluster_idx=1:numel(cluster_names)
    cluster_name=cluster_names(cluster_idx);
    
    %neighbors. from cluster "index" to neighbor "name"
    neigh{cluster_idx} = ...
        union(setdiff(cluster_annotation(2,cluster_annotation(1,:)==cluster_name),[0 cluster_name]), ...
              setdiff(cluster_annotation(1,cluster_annotation(2,:)==cluster_name),[0 cluster_name]));
    
    neigh_idx{cluster_idx} =  find(ismember(cluster_names,neigh{cluster_idx}));  
    
    %Indizen von nicht kopplenden Variablen im Cluster
    internal_var_indices{cluster_idx} = find( ...
        (cluster_id_from==cluster_name) .* (cluster_id_to==0) );
    
    %Indizen von nicht kopplenden + kopplenden Variablen im Cluster
    internal_plus_neighbor_var_indices{cluster_idx} = find( ...
        (cluster_id_from==cluster_name) + (cluster_id_to==cluster_name) );
    
    %Indizen von kopplenden Variablen im Cluster
    coupling_var_indices{cluster_idx} = ...
        intersect( internal_plus_neighbor_var_indices{cluster_idx}, ...
        find(cluster_id_to) );
    global_var_set = union(global_var_set, ...
        coupling_var_indices{cluster_idx});
    global_int_set = union(global_int_set, ...
        internal_var_indices{cluster_idx});
    
    relative_indices{cluster_idx} = [];
    for j=1:numel(coupling_var_indices{cluster_idx})
        relative_indices{cluster_idx} = [relative_indices{cluster_idx} ...
            find(global_var_set == coupling_var_indices{cluster_idx}(j))];
    end
    
    local_position_coupling_indices{cluster_idx} = [];
    for j=1:numel(coupling_var_indices{cluster_idx})
        local_position_coupling_indices{cluster_idx} = ...
            [local_position_coupling_indices{cluster_idx} ...
            find(internal_plus_neighbor_var_indices{cluster_idx} ...
            == coupling_var_indices{cluster_idx}(j))];
    end
    
    for neighbor=neigh{cluster_idx}
        %relative variable indices (from cluster "index" > neighbor "name")
        local_position_coupling_indices_with_neighbor{cluster_idx}{neighbor} = ...
            union(find(ismember(cluster_annotation(1,internal_plus_neighbor_var_indices{cluster_idx}),neighbor)),...
                  find(ismember(cluster_annotation(2,internal_plus_neighbor_var_indices{cluster_idx}),neighbor)));
    end
end

sub_models=cell(numel(cluster_names),1);
start_temp = cell(numel(cluster_names),1);

%Aufteilung A, lhs, rhs, ub, lb und f in Subprobleme
for cluster_idx=1:numel(cluster_names)
    
    [sub_models{cluster_idx}, f_sub{cluster_idx}, ...
        lin_cost_lambda{cluster_idx}, lin_cost_rho{cluster_idx}, ...
        diag_aux{cluster_idx}, constraint_rows{cluster_idx}] = ...
        initialize_subproblem(select, A, lhs, rhs, ub, lb, f, ... 
        cluster_idx, internal_plus_neighbor_var_indices, ...
        internal_var_indices, global_int_set); %#ok<AGROW>

end

%%Initializierung globalen Werte der kopplenden Variablen, lambda's und rho's
Plink_global = Plink_init * ones(column_size,1);
Plink_temp = zeros(column_size, numel(cluster_names));
ResidPrim_hist = zeros(1, maxiter);
ResidDual_hist = zeros(1, maxiter);
cost_hist = zeros(1, maxiter);
for cluster_idx = 1:numel(cluster_names)
    lambdas{cluster_idx} = lambdas_init * ...
        ones(size(coupling_var_indices{cluster_idx})); %#ok<AGROW>
    rhos{cluster_idx} = rhos_init; %#ok<AGROW>
    x_global{cluster_idx} = 0 * ones(size(coupling_var_indices{cluster_idx}));
    resid_prim{cluster_idx} = ...
        0 * ones(size(coupling_var_indices{cluster_idx}')); %#ok<AGROW>
    resid_dual{cluster_idx} = ...
        0 * ones(size(coupling_var_indices{cluster_idx}')); %#ok<AGROW>
    Plink_temp_coup{cluster_idx} = ...
        0*ones(size(coupling_var_indices{cluster_idx}')); %#ok<AGROW>
end

%Initializierung der ADMM-Klasse für asynkrone Jobs

for cluster_idx = 1:numel(cluster_names)
    sub{cluster_idx} = ADMMClass(sub_models{cluster_idx}, ...
        f_sub{cluster_idx}, ...
        lin_cost_lambda{cluster_idx}, ...
        lin_cost_rho{cluster_idx}, ...
        diag_aux{cluster_idx}, ...
        constraint_rows{cluster_idx},...
        lambdas{cluster_idx},...
        rhos{cluster_idx},...
        cluster_idx,...
        cluster_names,...
        neigh{cluster_idx},...
        internal_plus_neighbor_var_indices{cluster_idx}, ...
        internal_var_indices{cluster_idx}, ...
        local_position_coupling_indices{cluster_idx},...
        global_int_set,...
        neigh_idx{cluster_idx},...
        local_position_coupling_indices_with_neighbor{cluster_idx});
    sub{cluster_idx}.x_global = x_global{cluster_idx};  
end
solution=cell(numel(cluster_names));
%%ADMM-Iterationen
%Parallelpool erstellen (falls noch nicht vorhanden, oder vorhanden mit
%abweichender Threadanzahl)
try
    parpool(select.opt.num_threads)
catch
    parallel_cluster = parcluster;
    if parallel_cluster.NumWorkers ~= select.opt.num_threads
        delete(gcp('nocreate'));
        parpool(select.opt.num_threads)
    end
end


tic

spmd
    sub = g_runWorker(sub{labindex});
end


cpx = struct();

%write solution for main model
cpx.Solution.status = 1;
cpx.Solution.statusstring = 'optimal';
cpx.Solution.objval = 0;
cpx.Solution.dual = zeros(row_size,1);
cpx.Solution.ax = zeros(row_size,1);
cpx.Solution.reducedcost = zeros(column_size,1);
cpx.Solution.x = zeros(column_size,1);
cpx.Solution.time = toc;

for cluster_idx = 1:numel(cluster_names)
    sub_solution{cluster_idx} = sub{cluster_idx};
end
for cluster_idx = 1:numel(cluster_names)
    cpx.Solution.x(internal_plus_neighbor_var_indices{cluster_idx}) = sub_solution{cluster_idx}.solution;
end
cpx.Solution.objval = f*cpx.Solution.x;


display(strcat('Converged. ADMM.solution: ',num2str(cpx.Solution.x')))
display(strcat('Converged. Orig.solution: ',num2str(central_soln')))
display(strcat('Converged. ADMM.cost: ',num2str(cpx.Solution.objval)))
display(strcat('Converged. Orig.cost: ',num2str(central_soln_obj)))
end