function [cpx,central_soln] = solve_with_ADMM_parallel_Gabay_constant_rho_release(select, A, lhs, rhs, ub, lb, ...
    f, cluster_id_from, cluster_id_to)

cluster_id_from = cluster_id_from';
cluster_id_to = cluster_id_to';

row_size = size(A,1);
column_size = size(A,2);

%Build central model
orig_model = build_central(select, A, lhs, rhs, ub, lb, f);

%(optional) solve central problem to compare optimums
[central_soln, central_soln_obj] = linprog(orig_model.Model.obj, ...
                                           orig_model.Model.A, ...
                                           orig_model.Model.rhs, ...
                                           [], ...
                                           [], ...
                                           orig_model.Model.lb, ...
                                           orig_model.Model.ub);

%Read the ADMM parameters
maxiter = select.opt.decomposition.maxiter; 
mu = select.opt.decomposition.mu;
tau = select.opt.decomposition.tau;
Plink_init = select.opt.decomposition.Plink_init;
lambdas_init = select.opt.decomposition.lambdas_init;
rhos_init = select.opt.decomposition.rhos_init;
convergence_tolerance = select.opt.decomposition.convergence_tolerance;
convergence_type = select.opt.decomposition.convergence_type;

%Create an array with the cluster names
cluster_names = unique(cluster_id_from);

%Initialize and define index sets for each cluster
global_var_set = [];
global_int_set = [];
internal_var_indices = cell(numel(cluster_names),1);
internal_plus_neighbor_var_indices = cell(numel(cluster_names),1);
coupling_var_indices = cell(numel(cluster_names),1);
relative_indices = cell(numel(cluster_names),1);
local_position_coupling_indices = cell(numel(cluster_names),1);

for cluster_idx=1:numel(cluster_names)
    cluster_name=cluster_names(cluster_idx);
    
    %indices of non-coupling variables
    internal_var_indices{cluster_idx} = find( ...
        (cluster_id_from==cluster_name) .* (cluster_id_to==0) );
    
    %indices of non-coupling and coupling variables
    internal_plus_neighbor_var_indices{cluster_idx} = find( ...
        (cluster_id_from==cluster_name) + (cluster_id_to==cluster_name) );
    
    %indices of coupling variables
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
end

sub_models=cell(numel(cluster_names),1);
start_temp = cell(numel(cluster_names),1);

%Splitting of A, lhs, rhs, ub, lb und f into subproblems
for cluster_idx=1:numel(cluster_names)
    
    [sub_models{cluster_idx}, f_sub{cluster_idx}, ...
        lin_cost_lambda{cluster_idx}, lin_cost_rho{cluster_idx}, ...
        diag_aux{cluster_idx}, constraint_rows{cluster_idx}] = ...
        initialize_subproblem(select, A, lhs, rhs, ub, lb, f, ... 
        cluster_idx, internal_plus_neighbor_var_indices, ...
        internal_var_indices, global_int_set); %#ok<AGROW>
    
end

%Initizialization of some values
Plink_global = Plink_init * ones(column_size,1);
Plink_temp = zeros(column_size, numel(cluster_names));
ResidPrim_hist = zeros(1, maxiter);
ResidDual_hist = zeros(1, maxiter);
cost_hist = zeros(1, maxiter);

for cluster_idx = 1:numel(cluster_names)
    lambdas{cluster_idx} = lambdas_init * ...
        ones(size(coupling_var_indices{cluster_idx})); 
    rhos{cluster_idx} = rhos_init; 
    resid_prim{cluster_idx} = ...
        0 * ones(size(coupling_var_indices{cluster_idx}'));
    resid_dual{cluster_idx} = ...
        0 * ones(size(coupling_var_indices{cluster_idx}')); 
    Plink_temp_coup{cluster_idx} = ...
        0*ones(size(coupling_var_indices{cluster_idx}')); 
end

solution=cell(numel(cluster_names),1);
local_cost=cell(numel(cluster_names),1);


%create a parallel pool
try
    parpool(select.opt.num_threads)
catch
    parallel_cluster = parcluster;
    if parallel_cluster.NumWorkers ~= select.opt.num_threads
        delete(gcp('nocreate'));
        parpool(select.opt.num_threads)
    end
end

%Figures: 1) Cost  function, 3) Primal residual, 4) Dual residual 
kosten_vergleich=subplot(221)
primal_residuals_progress=subplot(223)
dual_residuals_progress=subplot(224)

tic
%%ADM Iterations
for iter=1:maxiter
    display(['Start iteration ',num2str(iter)]);
    
    
    %Relevante globale Werte für alle Cluster zuweisen
    for cluster_idx=1:numel(cluster_names)
        Plink_global_parallel{cluster_idx} = ...
            Plink_global(coupling_var_indices{cluster_idx})'; %#ok<AGROW>
        x{cluster_idx} = 0 * ones(size(internal_plus_neighbor_var_indices{cluster_idx}'));
    end    
    
    %Parallel-Loop
    parfor cluster_idx=1:numel(cluster_names)
        disp(['Building cluster ',num2str(cluster_names(cluster_idx)), ...
            ' (',num2str(cluster_idx),'th cluster)']);
        
        %linear penalty terms from lambda und rho
        lin_cost_lambda{cluster_idx}(local_position_coupling_indices{cluster_idx}) = lambdas{cluster_idx}; %#ok<PFOUS>
        lin_cost_rho{cluster_idx}(local_position_coupling_indices{cluster_idx}) =  -rhos{cluster_idx} .* Plink_global_parallel{cluster_idx}; %#ok<PFOUS>
        
        %linear penalty terms are added to the original cost function
        f_sub_augmented{cluster_idx} = f_sub{cluster_idx} + ...
            lin_cost_lambda{cluster_idx} + lin_cost_rho{cluster_idx}; %#ok<PFOUS>
        
        %quadratic penalty terms from rho
        diag_aux{cluster_idx}(local_position_coupling_indices{cluster_idx}) ...
            = rhos{cluster_idx}; %#ok<PFOUS>
        h_sub{cluster_idx}=spdiags(diag_aux{cluster_idx}, 0, ...
            size(internal_plus_neighbor_var_indices{cluster_idx}, 1), ...
            size(internal_plus_neighbor_var_indices{cluster_idx}, 1)); %#ok<PFOUS>
        
        sub_models{cluster_idx}.Model.obj = f_sub_augmented{cluster_idx}; %#ok<PFOUS>
        sub_models{cluster_idx}.Model.Q = h_sub{cluster_idx};

        %solve subproblem
        disp(['Iteration #',num2str(iter),': Solving cluster ', ...
            num2str(cluster_names(cluster_idx)),' (', ...
            num2str(cluster_idx),'th cluster)']);
                              
        [solution{cluster_idx}, local_cost{cluster_idx}] = quadprog(sub_models{cluster_idx}.Model.Q,...
                                         sub_models{cluster_idx}.Model.obj,...
                                         sub_models{cluster_idx}.Model.A,...
                                         sub_models{cluster_idx}.Model.rhs,...
                                         [],[],...
                                         sub_models{cluster_idx}.Model.lb,...
                                         sub_models{cluster_idx}.Model.ub);
                                     

        disp(['Iteration #',num2str(iter),': Solved cluster ', ...
            num2str(cluster_names(cluster_idx)),' (', ...
            num2str(cluster_idx),'th cluster)']);
        
        x{cluster_idx} = solution{cluster_idx};
        
        %difference between two iterations calculated for the dual residual
        if iter>1
            x_diff_squared{cluster_idx} = ...
                (x{cluster_idx}(local_position_coupling_indices{cluster_idx}) - ...
                Plink_temp_coup{cluster_idx}).^2;
        end
        
        Plink_temp_parallel{cluster_idx} = x{cluster_idx};
    end
    
    for cluster_idx=1:numel(cluster_names)
        Plink_temp(internal_plus_neighbor_var_indices{cluster_idx}, ...
            cluster_idx) = Plink_temp_parallel{cluster_idx};
        Plink_temp_coup{cluster_idx} = ...
            Plink_temp(coupling_var_indices{cluster_idx}, cluster_idx);
    end
    
    %average the global values
    P_link_summed = sum(Plink_temp,2);
    P_link_summed(global_var_set) = P_link_summed(global_var_set) / 2;
    Plink_global = P_link_summed;
    cost_hist(iter) = f*Plink_global;
    
    if iter>1
        for cluster_idx=1:numel(cluster_names)
            resid_prim{cluster_idx} = norm( ...
                (x{cluster_idx}(local_position_coupling_indices{cluster_idx}) - ...
                Plink_global(coupling_var_indices{cluster_idx})));
            resid_dual{cluster_idx} = rhos{cluster_idx}(1) * ...
                sqrt(sum(x_diff_squared{cluster_idx}));
            
        end
        
    end
    
    %update lambda
    for cluster_idx = 1:numel(cluster_names)
        lambdas{cluster_idx} = lambdas{cluster_idx} + rhos{cluster_idx} .* ...
            (x{cluster_idx}(local_position_coupling_indices{cluster_idx}) - ...
            Plink_global(coupling_var_indices{cluster_idx}));
    end
        
    %check convergence
    if iter > 1
        ResidPrim_hist(iter) = norm(cell2mat(resid_prim));
        ResidDual_hist(iter) = norm(cell2mat(resid_dual));
        hold on
        if convergence_type == 1
            if norm(cell2mat(resid_prim)) < convergence_tolerance && ...
                    norm(cell2mat(resid_dual)) < convergence_tolerance
                break
            end
        elseif convergence_type == 2
            if abs(central_soln_obj - f*Plink_global)/abs(central_soln_obj) < convergence_tolerance
                break
            end
        end
    end
    disp(['Iteration #',num2str(iter), '. ADMM.cost: ', ...
        num2str(f*Plink_global), '. Iteration duration: ', ...
        num2str(toc), ' seconds.'])
end

%write solution for main model
cpx.Solution.status = 1;
cpx.Solution.statusstring = 'optimal';
cpx.Solution.objval = 0;
cpx.Solution.x = zeros(column_size,1);
cpx.Solution.time = toc;

for cluster_idx = 1:numel(cluster_names)
    cpx.Solution.objval = cpx.Solution.objval + local_cost{cluster_idx};
    cpx.Solution.x(internal_plus_neighbor_var_indices{cluster_idx}) = solution{cluster_idx};
end
semilogy(primal_residuals_progress,ResidPrim_hist(1:iter)')
semilogy(dual_residuals_progress,ResidDual_hist(1:iter)')
plot(kosten_vergleich,1:iter,cost_hist(1:iter),'x',1:iter,central_soln_obj*ones(size(1:iter)),'o')
display(strcat('Converged. ADMM.solution: ',num2str(Plink_global')))
display(strcat('Converged. Orig.solution: ',num2str(central_soln')))
display(strcat('Converged. ADMM.cost: ',num2str(f*Plink_global)))
display(strcat('Converged. Orig.cost: ',num2str(f*central_soln)))
end