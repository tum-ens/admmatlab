function [cpx, central_soln] = solve_with_ADMM_sequential_Song_local_rho_release(select, A, lhs, rhs, ub, lb, ...
    f, cluster_id_from, cluster_id_to)

cluster_id_from = cluster_id_from';
cluster_id_to = cluster_id_to';

row_size = size(A,1);
column_size = size(A,2);

%Zentrales Problem aufbauen f�r die Output-Struct (cpx)
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
iter_rho_adapt = select.opt.decomposition.iter_rho_adapt;
convergence_type = select.opt.decomposition.convergence_type;

%tau = 0;
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
%Plink_hist = zeros(size(A,2), maxiter);
ResidPrim_hist = zeros(1, maxiter);
ResidDual_hist = zeros(1, maxiter);
cost_hist = zeros(1, maxiter);

for cluster_idx = 1:numel(cluster_names)
    lambdas{cluster_idx} = lambdas_init * ...
        ones(size(coupling_var_indices{cluster_idx})); %#ok<AGROW>
    rhos{cluster_idx} = rhos_init; %#ok<AGROW>
    resid_prim{cluster_idx} = ...
        0 * ones(size(coupling_var_indices{cluster_idx}')); %#ok<AGROW>
    resid_dual{cluster_idx} = ...
        0 * ones(size(coupling_var_indices{cluster_idx}')); %#ok<AGROW>
    Plink_temp_coup{cluster_idx} = ...
        0*ones(size(coupling_var_indices{cluster_idx}')); %#ok<AGROW>
end

solution=cell(numel(cluster_names),1);
local_cost=cell(numel(cluster_names),1);

%%ADMM-Iterationen
%Zwei Figuren: 1) Kostenfunktion, 3) Primal residual, 4) Dual residual
kosten_vergleich=subplot(221)
primal_residuals_progress=subplot(223)
dual_residuals_progress=subplot(224)

tic
for iter=1:maxiter
    display(['Start iteration ',num2str(iter)]);
    
    %Relevante globale Werte f�r alle Cluster zuweisen
    for cluster_idx=1:numel(cluster_names)
        Plink_global_parallel{cluster_idx} = ...
            Plink_global(coupling_var_indices{cluster_idx})'; %#ok<AGROW>
        x{cluster_idx} = 0 * ones(size(internal_plus_neighbor_var_indices{cluster_idx}'));
    end
    
    %Parallel-Loop
    for cluster_idx=1:numel(cluster_names)
        disp(['Building cluster ',num2str(cluster_names(cluster_idx)), ...
            ' (',num2str(cluster_idx),'th cluster)']);
        
        %Lineare Kostenterme aus lambda und rho
        lin_cost_lambda{cluster_idx}(local_position_coupling_indices{cluster_idx}) = lambdas{cluster_idx}'; %#ok<PFOUS>
        lin_cost_rho{cluster_idx}(local_position_coupling_indices{cluster_idx}) =  -rhos{cluster_idx} .* Plink_global_parallel{cluster_idx}; %#ok<PFOUS>
        
        %lineare Kostenterme addiert auf die urspr�ngliche Kostenfunktion
        f_sub_augmented{cluster_idx} = f_sub{cluster_idx} + ...
            lin_cost_lambda{cluster_idx} + lin_cost_rho{cluster_idx}; %#ok<PFOUS>
        
        %quadratischer Kostenterm aus rho
        diag_aux{cluster_idx}(local_position_coupling_indices{cluster_idx}) ...
            = rhos{cluster_idx}; %#ok<PFOUS>
        h_sub{cluster_idx}=spdiags(diag_aux{cluster_idx}, 0, ...
            size(internal_plus_neighbor_var_indices{cluster_idx}, 1), ...
            size(internal_plus_neighbor_var_indices{cluster_idx}, 1)); %#ok<PFOUS>
        
        sub_models{cluster_idx}.Model.obj = f_sub_augmented{cluster_idx}; %#ok<PFOUS>
        sub_models{cluster_idx}.Model.Q = h_sub{cluster_idx};
        
        %Subprobleme l�sen
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
        
        %Differenz zwischen Iterationen berechnet f�r Dualresidual
        if iter>1
            x_diff{cluster_idx} = ...
                (x{cluster_idx}(local_position_coupling_indices{cluster_idx}) - ...
                Plink_temp_coup{cluster_idx});
        end
        
        Plink_temp_parallel{cluster_idx} = x{cluster_idx};
    end
    
    for cluster_idx=1:numel(cluster_names)
        Plink_temp(internal_plus_neighbor_var_indices{cluster_idx}, ...
            cluster_idx) = Plink_temp_parallel{cluster_idx};
        Plink_temp_coup{cluster_idx} = ...
            Plink_temp(coupling_var_indices{cluster_idx}, cluster_idx);
    end
    
    %Globale Werte mitteln
    P_link_summed = sum(Plink_temp,2);
    P_link_summed(global_var_set) = P_link_summed(global_var_set) / 2;
    Plink_global = P_link_summed;
    cost_hist(iter) = f*Plink_global;
    
    % penalty parameter tuning
    if iter>1
        for cluster_idx=1:numel(cluster_names)
            resid_prim{cluster_idx} = x{cluster_idx}(local_position_coupling_indices{cluster_idx}) - ...
                Plink_global(coupling_var_indices{cluster_idx});
            resid_dual{cluster_idx} = rhos{cluster_idx} * ...
                x_diff{cluster_idx};
            
            
            
            % cluster specific penalty parameters (no convergence
            % guarantee, but makes faster)
            idx_rho_increase = (norm(resid_prim{cluster_idx}) > ...
                mu * norm(resid_dual{cluster_idx}));
            idx_rho_decrease = (norm(resid_dual{cluster_idx}) > ...
                mu * norm(resid_prim{cluster_idx}));
            idx_rho_equal = ~(idx_rho_increase | idx_rho_decrease);
            
            if iter < iter_rho_adapt
                rhos{cluster_idx} = idx_rho_increase * ...
                    rhos{cluster_idx} * (1 + tau) + ...
                    idx_rho_decrease * rhos{cluster_idx} / (1 + tau) + ...
                    idx_rho_equal * rhos{cluster_idx};
                
            end
        end
        
        %Lambda updaten
        for cluster_idx = 1:numel(cluster_names)
            lambdas{cluster_idx} = lambdas{cluster_idx} + rhos{cluster_idx} .* ...
                (x{cluster_idx}(local_position_coupling_indices{cluster_idx}) - ...
                Plink_global(coupling_var_indices{cluster_idx}));
        end
        
        plot(kosten_vergleich,1:iter,cost_hist(1:iter),'x',1:iter,central_soln_obj*ones(size(1:iter)),'o')
        ylim(kosten_vergleich,[(min(cost_hist(iter),central_soln_obj)-abs(cost_hist(iter)-central_soln_obj)),...
            (max(cost_hist(iter),central_soln_obj)+abs(cost_hist(iter)-central_soln_obj))])
        
        pause(0)
        
        %Konvergenz pr�fen
        if iter > 1
            ResidPrim_hist_sum(iter) = 0;
            ResidDual_hist_sum(iter) = 0;
            for cluster_idx = 1:numel(cluster_names)
                ResidPrim_hist_sum(iter) = ResidPrim_hist_sum(iter) + norm(resid_prim{cluster_idx}) ;
                ResidDual_hist_sum(iter) = ResidDual_hist_sum(iter) + norm(resid_dual{cluster_idx});
            end
            semilogy(primal_residuals_progress,ResidPrim_hist_sum(1:iter)')
            semilogy(dual_residuals_progress,ResidDual_hist_sum(1:iter)')
            hold on
            if convergence_type == 1
                if ResidPrim_hist_sum(iter) < convergence_tolerance && ...
                        ResidDual_hist_sum(iter) < convergence_tolerance
                    break
                end
            elseif convergence_type == 2
                if abs(central_soln_obj - f*Plink_global)/abs(central_soln.objval) < convergence_tolerance
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
    cpx.Solution.dual = zeros(row_size,1);
    cpx.Solution.ax = zeros(row_size,1);
    cpx.Solution.reducedcost = zeros(column_size,1);
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