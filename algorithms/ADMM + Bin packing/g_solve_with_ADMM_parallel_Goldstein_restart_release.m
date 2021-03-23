function [cpx,central_soln] = g_solve_with_ADMM_parallel_Goldstein_restart_release(select, A, lhs, rhs, ub, lb, ...
    f, cluster_id_from, cluster_id_to)

cluster_id_from = cluster_id_from';
cluster_id_to = cluster_id_to';

row_size = size(A,1);
column_size = size(A,2);

%Zentrales Problem aufbauen
orig_model = build_central(select, A, lhs, rhs, ub, lb, f);

%Bin Packing
bin_size = select.opt.bin_size; % how many x8 per bin, perhaps via input
bin_number = floor(select.opt.num_threads/bin_size); % number of bins

%define sub_cluster_masks
for i = 1:bin_number-1
    subnumthreads{i} = bin_size;
end

subnumthreads{bin_number} = select.opt.num_threads-(bin_size*(bin_number-1));

%parallel setting fuer Vergleichsproblem
orig_model.Param.threads.Cur = subnumthreads{1};

if isunix
    orig_model.Param.cpumask.Cur = submask{1}.cplex;
end
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
mu_sc = select.opt.decomposition.convergence_tolerance;
convergence_type = select.opt.decomposition.convergence_type;

alpha = ones(1,maxiter+1);
c = zeros(1,maxiter+1);
c(1) = 1;

%Ein Array von eindeutigen Clusternamen
cluster_names = unique(cluster_id_from);

%Verschiedene Indizmengen definieren f�r jeweilige Cluster
global_var_set = [];
global_int_set = [];
internal_var_indices = cell(numel(cluster_names),1);
internal_plus_neighbor_var_indices = cell(numel(cluster_names),1);
coupling_var_indices = cell(numel(cluster_names),1);
relative_indices = cell(numel(cluster_names),1);
local_position_coupling_indices = cell(numel(cluster_names),1);

for cluster_idx=1:numel(cluster_names)
    cluster_name=cluster_names(cluster_idx);
    
    %Indizen von nicht kopplenden Variablen im Cluster
    internal_var_indices{cluster_idx} = find( (cluster_id_from==cluster_name) .* (cluster_id_to==0) );
    
    %Indizen von nicht kopplenden + kopplenden Variablen im Cluster
    internal_plus_neighbor_var_indices{cluster_idx} = find( (cluster_id_from==cluster_name) + (cluster_id_to==cluster_name) );
    
    %Indizen von kopplenden Variablen im Cluster
    coupling_var_indices{cluster_idx} = intersect( internal_plus_neighbor_var_indices{cluster_idx}, find(cluster_id_to) );
    global_var_set = union(global_var_set, coupling_var_indices{cluster_idx});
    global_int_set = union(global_int_set, internal_var_indices{cluster_idx});
    
    relative_indices{cluster_idx} = [];
    for j=1:numel(coupling_var_indices{cluster_idx})
        relative_indices{cluster_idx} = [relative_indices{cluster_idx} find(global_var_set == coupling_var_indices{cluster_idx}(j))];
    end
    
    local_position_coupling_indices{cluster_idx} = [];
    for j=1:numel(coupling_var_indices{cluster_idx})
        local_position_coupling_indices{cluster_idx} = [local_position_coupling_indices{cluster_idx} find(internal_plus_neighbor_var_indices{cluster_idx} == coupling_var_indices{cluster_idx}(j))];
    end
end

sub_models=cell(numel(cluster_names),1);
start_temp = cell(numel(cluster_names),1);

%Aufteilung A, lhs, rhs, ub, lb und f in Subprobleme
for cluster_idx=1:numel(cluster_names)
    cluster_name=cluster_names(cluster_idx);
    
    [sub_models{cluster_idx}, f_sub{cluster_idx}, ...
        lin_cost_lambda{cluster_idx}, lin_cost_rho{cluster_idx}, ...
        diag_aux{cluster_idx}, constraint_rows{cluster_idx}] = ...
        initialize_subproblem(select, A, lhs, rhs, ub, lb, f, ...
        cluster_idx, internal_plus_neighbor_var_indices, ...
        internal_var_indices, global_int_set); %#ok<AGROW>
    
end

% weight for Bin Packing
weights = zeros(1,numel(cluster_names));
for cluster_idx = 1:numel(cluster_names)
    weights(cluster_idx) = size(sub_models{cluster_idx}.Model.A,2);
end

% Bin Packing via LPT (Longest Processing Time)
grouping = group_subproblems_alternative(weights,bin_number);
% masks per bin
for bin = 1:bin_number
    for idx = grouping{bin}
        sub_models{idx}.Param.threads.Cur = subnumthreads{bin};
        if isunix
            sub_models{idx}.Param.cpumask.Cur = submask{bin}.cplex;
        end
    end
end

%%Initializierung globalen Werte der kopplenden Variablen, lambda's und rho's
Plink_global = Plink_init * ones(column_size,maxiter+1);
Plink_global_hat = Plink_init * ones(column_size,maxiter);

Plink_temp = zeros(column_size, numel(cluster_names));
ResidPrim_hist = zeros(1, maxiter);
ResidDual_hist = zeros(1, maxiter);

for cluster_idx = 1:numel(cluster_names)
    lambdas{cluster_idx} = lambdas_init * ones(size(coupling_var_indices{cluster_idx},1),maxiter+1);
    lambdas_hat{cluster_idx} = lambdas_init * ones(size(coupling_var_indices{cluster_idx},1),maxiter);
    rhos{cluster_idx} = rhos_init; %#ok<AGROW>
    resid_prim{cluster_idx} = ...
        0 * ones(size(coupling_var_indices{cluster_idx}')); %#ok<AGROW>
    resid_dual{cluster_idx} = ...
        0 * ones(size(coupling_var_indices{cluster_idx}')); %#ok<AGROW>
    Plink_temp_coup{cluster_idx} = ...
        0*ones(size(coupling_var_indices{cluster_idx}')); %#ok<AGROW>end
end
solution_in_bins=cell(1,bin_number);
solution=cell(1,numel(cluster_names));

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
kosten_vergleich=subplot(221)
primal_residuals_progress=subplot(223)
dual_residuals_progress=subplot(224)
set(dual_residuals_progress,'yscale','log')
set(primal_residuals_progress,'yscale','log')

tic
for iter=1:maxiter
    time_until_now = toc;
    display(['Start iteration ',num2str(iter)]);
    %Relevante globale Werte f�r alle Cluster zuweisen
    for cluster_idx=1:numel(cluster_names)
        Plink_global_parallel{cluster_idx} = ...
            Plink_global_hat(coupling_var_indices{cluster_idx},iter)'; %#ok<AGROW>
        if iter==1
            x{cluster_idx} = 0 * ones(size(internal_plus_neighbor_var_indices{cluster_idx}'));
        end
        lambdas_hat_curr{cluster_idx}=lambdas_hat{cluster_idx}(:,iter);
    end
    
    for cluster_idx=1:numel(cluster_names)
        %disp(['Building cluster ',num2str(cluster_names(cluster_idx)),' (',num2str(cluster_idx),'th cluster)']);
        lin_cost_lambda{cluster_idx}(local_position_coupling_indices{cluster_idx}) = lambdas_hat_curr{cluster_idx};
        lin_cost_rho{cluster_idx}(local_position_coupling_indices{cluster_idx}) =  -rhos{cluster_idx} .* Plink_global_parallel{cluster_idx};
        
        %lineare Kostenterme addiert auf die urspr�ngliche Kostenfunktion
        sub_models{cluster_idx}.Model.obj = f_sub{cluster_idx} + lin_cost_lambda{cluster_idx} + lin_cost_rho{cluster_idx};
        
        %sub_models{cluster_idx}.Model.obj = f_sub_augmented{cluster_idx};
        
        %quadratischer Kostenterm aus rho
        diag_aux{cluster_idx}(local_position_coupling_indices{cluster_idx}) = rhos{cluster_idx};
        sub_models{cluster_idx}.Model.Q = spdiags(diag_aux{cluster_idx}, 0,...
            size(internal_plus_neighbor_var_indices{cluster_idx},1), ...
            size(internal_plus_neighbor_var_indices{cluster_idx},1));
        
    end
    %solve subproblems as bins in parallel
    for bin=1:bin_number
        for cluster_idx=grouping{bin}
            sub_models_grouped{bin,cluster_idx}=sub_models{cluster_idx};
        end
    end
    parfor bin=1:bin_number
        solution_temp = cell(1,numel(grouping{bin}));
        for cluster_idx=grouping{bin}
            model_temp = sub_models_grouped{bin,cluster_idx};
            
            solution_temp{cluster_idx} = quadprog(model_temp.Model.Q,...
                model_temp.Model.obj,...
                model_temp.Model.A,...
                model_temp.Model.rhs,...
                [],[],...
                model_temp.Model.lb,...
                model_temp.Model.ub);
        end
        solution_in_bins{bin} = solution_temp;
    end
    
    %recover subproblems
    for bin=1:bin_number
        for cluster_idx=grouping{bin}
            solution{cluster_idx}=solution_in_bins{bin}{cluster_idx};
        end
    end
    for cluster_idx=1:numel(cluster_names)
        x{cluster_idx} = solution{cluster_idx};
        
        if iter>1
            x_diff_squared{cluster_idx} = (x{cluster_idx}(local_position_coupling_indices{cluster_idx}) - ...
                Plink_temp_coup{cluster_idx}).^2;
        end
        Plink_temp_parallel{cluster_idx} = x{cluster_idx};
    end
    for cluster_idx=1:numel(cluster_names)
        %         sub_models{cluster_idx}.Start = start_temp{cluster_idx};
        Plink_temp(internal_plus_neighbor_var_indices{cluster_idx}, ...
            cluster_idx) = Plink_temp_parallel{cluster_idx};
        Plink_temp_coup{cluster_idx} = ...
            Plink_temp(coupling_var_indices{cluster_idx}, cluster_idx);
    end
    
    %Globale Werte mitteln
    P_link_summed = sum(Plink_temp,2);
    P_link_summed(global_var_set) = P_link_summed(global_var_set) / 2;
    Plink_global(:,iter+1) = P_link_summed;
    Plink_hist(:, iter) = Plink_global(:,iter+1);
    cost_hist(iter) = f*Plink_global(:,iter+1);
    
    % penalty parameter tuning
    if iter>1
        for cluster_idx=1:numel(cluster_names)
            resid_prim{cluster_idx}=norm((x{cluster_idx}(local_position_coupling_indices{cluster_idx}) - ...
                Plink_global(coupling_var_indices{cluster_idx},iter+1)));
            resid_dual{cluster_idx}=rhos{cluster_idx}(1)*sqrt(sum(x_diff_squared{cluster_idx}));
            
        end
        % same penalty parameter for all clusters
        idx_rho_increase = (norm(cell2mat(resid_prim)) > mu * norm(cell2mat(resid_dual)));
        idx_rho_decrease = (norm(cell2mat(resid_prim)) < mu * norm(cell2mat(resid_dual)));
        idx_rho_equal = (norm(cell2mat(resid_prim)) == mu * norm(cell2mat(resid_dual)));
        
        
        for cluster_idx=1:numel(cluster_names)
            
            rhos{cluster_idx}=idx_rho_increase .* rhos{cluster_idx}*(1+tau)+ ...
                idx_rho_decrease .* rhos{cluster_idx}/(1+tau)+ ...
                idx_rho_equal .* rhos{cluster_idx};
        end
    end
    
    %Lambda updaten
    for cluster_idx = 1:numel(cluster_names)
        lambdas{cluster_idx}(:,iter+1) = lambdas_hat{cluster_idx}(:,iter) + rhos{cluster_idx} .* ...
            (x{cluster_idx}(local_position_coupling_indices{cluster_idx}) - Plink_global(coupling_var_indices{cluster_idx},iter+1));
    end
    
    c(iter+1)=0;
    %Calculating the combined residual
    for cluster_idx = 1:numel(cluster_names)
        c(iter+1) = c(iter+1)+ ...
            rhos{1}(1)^-1 / (numel(cluster_names))^2 * norm(lambdas{cluster_idx}(:,iter+1)-lambdas_hat{cluster_idx}(:,iter))^2;
    end
    
    c(iter+1) = c(iter+1) + rhos{1}(1) * norm(Plink_global(:,iter+1)-Plink_global_hat(:,iter))^2;
    
    if c(iter+1) < mu_sc * c(iter)
        alpha(iter+1) = (1+sqrt(1+4*alpha(iter)^2))/2;
        %to disable the fast ADMM
        %alpha(iter+1) = 1;
        Plink_global_hat(:,iter+1) = Plink_global(:,iter+1) + (alpha(iter)-1)/(alpha(iter+1))*(Plink_global(:,iter+1) - Plink_global(:,iter));
        
        for cluster_idx = 1:numel(cluster_names)
            lambdas_hat{cluster_idx}(:,iter+1) = lambdas{cluster_idx}(:,iter+1) + (alpha(iter)-1)/(alpha(iter+1))*(lambdas{cluster_idx}(:,iter+1)-lambdas{cluster_idx}(:,iter));
        end
    else
        alpha(iter+1) = 1;
        Plink_global_hat(:,iter+1) = Plink_global(:,iter);
        for cluster_idx = 1:numel(cluster_names)
            lambdas_hat{cluster_idx}(:,iter+1) = lambdas{cluster_idx}(:,iter);
        end
    end
    
    plot(kosten_vergleich,1:iter,cost_hist(1:iter),'x',1:iter,central_soln_obj*ones(size(1:iter)),'o')
    ylim(kosten_vergleich,sort([min(cost_hist(iter),central_soln_obj)*0.9,max(cost_hist(iter),central_soln_obj)*1.1]))
    pause(0)
    if iter > 1
        ResidPrim_hist(iter) = norm(cell2mat(resid_prim));
        ResidDual_hist(iter) = norm(cell2mat(resid_dual));
        semilogy(primal_residuals_progress,1:iter,ResidPrim_hist(1:iter)')
        
        semilogy(dual_residuals_progress,1:iter,ResidDual_hist(1:iter)')
        
        if convergence_type == 1
            if norm(cell2mat(resid_prim)) < convergence_tolerance && ...
                    norm(cell2mat(resid_dual)) < convergence_tolerance
                break
            end
        elseif convergence_type == 2
            if abs(central_soln.objval - f*num2str(f*Plink_global(:,iter)))/central_soln_obj < convergence_tolerance
                break
            end
        end
    end
    disp(['Iteration #',num2str(iter), '. ADMM.cost: ',num2str(f*Plink_global(:,iter)), ...
        ' Iteration took ',num2str(toc-time_until_now),' seconds'])
end

cpx.Solution.status = 1;
cpx.Solution.statusstring = 'optimal';
cpx.Solution.objval = 0;
cpx.Solution.dual = zeros(row_size,1);
cpx.Solution.ax = zeros(row_size,1);
cpx.Solution.reducedcost = zeros(column_size,1);
cpx.Solution.x = zeros(column_size,1);
cpx.Solution.time = toc;


for cluster_idx = 1:numel(cluster_names)
    cpx.Solution.x(internal_plus_neighbor_var_indices{cluster_idx}) = solution{cluster_idx};
end
cpx.iterations = iter;
cpx.Solution.objval = f*Plink_global;

semilogy(primal_residuals_progress,ResidPrim_hist(1:iter)')
semilogy(dual_residuals_progress,ResidDual_hist(1:iter)')
plot(kosten_vergleich,1:iter,cost_hist(1:iter),'x',1:iter,central_soln_obj*ones(size(1:iter)),'o')
display(strcat('Converged. ADMM.solution: ',num2str(Plink_global(:,iter)')))
display(strcat('Converged. Orig.solution: ',num2str(central_soln')))
display(strcat('Converged. ADMM.cost: ',num2str(f*Plink_global(:,iter))))
display(strcat('Converged. Orig.cost: ',num2str(f*central_soln)))

end