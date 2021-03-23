function obj = g_runWorker(obj)
maxit = 1000;
resid_prim = 999;

tol = 10^-5;

rho_max = 10;
tau = 1.1;
tau_max = 100;
mu = 1;
obj.rhos = 0.1;
eta = 1;

tag_from_neighbors=ones([1, numel(obj.neigh_idx)]);
most_recent_data_from_neighbors = zeros([1, numel(obj.neigh_idx)]);
convergences = zeros([1, numel(obj.neigh_idx)]);

rho_from_neighbors = zeros([1, numel(obj.neigh_idx)]);

resid_prim_with_neighbor=zeros([1, numel(obj.neigh_idx)]);
resid_dual_with_neighbor=zeros([1, numel(obj.neigh_idx)]);
converge = 0;

convergence_allowance_to_neighbor = zeros([1, numel(obj.neigh_idx)]);

for iter = 1:maxit
    new_msg_received=zeros([1, numel(obj.neigh_idx)]);
    convergence_allowance_from_neighbor = zeros([1, numel(obj.neigh_idx)]);
    
    %% if recvmsg: update f_global for line connecting
    if iter>1
        %% dont act until a neighbor finished solving
        if numel(setdiff(obj.neigh, obj.neigh(find(convergences))))>0
            while sum(new_msg_received) < 1
                for neighbor_name = setdiff(obj.neigh, obj.neigh(find(convergences)))
                    if labProbe(find(obj.cluster_names==neighbor_name), maxit*100*find(obj.cluster_names==neighbor_name))
                        tag_for_data(obj.neigh==neighbor_name) = labReceive(find(obj.cluster_names==neighbor_name), maxit*100*find(obj.cluster_names==neighbor_name));
                        [new_msg_received(obj.neigh==neighbor_name), ~, tag_]  = ...
                            labProbe(find(obj.cluster_names==neighbor_name), tag_for_data(obj.neigh==neighbor_name));
                    end
                end
            end
        end
        %% check convergence and get rhos from neighbors
        for rel_neig_index = find(convergences)
            convergence_allowance_from_neighbor(rel_neig_index)= 1;
        end
        
        for rel_neig_index = find(new_msg_received)
            neighbor = obj.neigh(rel_neig_index);
            convergences(rel_neig_index) = labReceive(find(obj.cluster_names==obj.neigh(rel_neig_index)), ...
                maxit*200*find(obj.cluster_names==obj.neigh(rel_neig_index)));
            convergence_allowance_from_neighbor(rel_neig_index) = labReceive(find(obj.cluster_names==obj.neigh(rel_neig_index)), ...
                maxit*400*find(obj.cluster_names==obj.neigh(rel_neig_index)));
            convergence_allowance_from_neighbor(rel_neig_index)= 1;
            
            if rem(iter,10)==0
                disp(['Iteration #',num2str(iter),': Cluster ', ...
                    num2str(obj.cluster_names(labindex)), ' . Neighbors: ',...
                    num2str(obj.neigh), ...
                    'convergence allowances from neighbors: ',num2str(convergence_allowance_from_neighbor),...
                    'Self convergence?: ', num2str(converge)])
            end
        end
        
        if converge
            disp(['Cluster with index ',num2str(labindex), ' and name ', num2str(obj.cluster_names(labindex)), ...
                ' converged at the ',num2str(iter) ,'. iteration. ', ...
                'Neighbor convergences: ',num2str(convergences),...
                'Self convergence?: ', num2str(converge)]);
            obj.it = iter;
            break
        end
        for rel_neig_index=find(new_msg_received)
            data_from_neighbor{rel_neig_index} = labReceive(find(obj.cluster_names==obj.neigh(rel_neig_index)), tag_for_data(rel_neig_index));
            
            tag_from_neighbors(rel_neig_index) = tag_from_neighbors(rel_neig_index) + 1;
            neighbor = obj.neigh(rel_neig_index);
            
            f_global_hist = obj.x_global;
            obj.x_global(...
                ismember(obj.local_position_coupling_indices,...
                obj.local_position_coupling_indices_with_neighbor{neighbor})) = ...
                (obj.solution(obj.local_position_coupling_indices_with_neighbor{neighbor})'+ ...
                data_from_neighbor{rel_neig_index}')/2;
        end
    end
    
    %Lineare Kostenterme aus lambda und rho
    obj.lin_cost_lambda(obj.local_position_coupling_indices) = obj.lambdas; %#ok<PFOUS>
    obj.lin_cost_rho(obj.local_position_coupling_indices) =  -obj.rhos * obj.x_global; %#ok<PFOUS>
    
    %lineare Kostenterme addiert auf die urspr�ngliche Kostenfunktion
    f_sub_augmented = obj.f_sub + ...
        obj.lin_cost_lambda + obj.lin_cost_rho; %#ok<PFOUS>
    
    %quadratischer Kostenterm aus rho
    obj.diag_aux(obj.local_position_coupling_indices) ...
        = obj.rhos; %#ok<PFOUS>
    h_sub=spdiags(obj.diag_aux, 0, ...
        size(obj.internal_plus_neighbor_var_indices, 1), ...
        size(obj.internal_plus_neighbor_var_indices, 1)); %#ok<PFOUS>
    
    obj.sub_model.Model.obj = f_sub_augmented; %#ok<PFOUS>
    obj.sub_model.Model.Q = h_sub;
    
    %Subprobleme l�sen
    obj.solution = quadprog(obj.sub_model.Model.Q,...
        obj.sub_model.Model.obj,...
        obj.sub_model.Model.A,...
        obj.sub_model.Model.rhs,...
        [],[],...
        obj.sub_model.Model.lb,...
        obj.sub_model.Model.ub);...,X0,options,varargin)
        
    x = obj.solution;
    
    %% update lambda
    if iter>1
        % update all lambdas
        obj.lambdas = ...
            obj.lambdas +...
            obj.rhos *  (x(obj.local_position_coupling_indices) - ...
            obj.x_global);
    end
    
    %% update rho
    if iter>1
        %residuals with neighbor
        for rel_neig_index=find(new_msg_received)
            %relative residuals
            neighbor = obj.neigh(rel_neig_index);
            
            resid_prim_with_neighbor(rel_neig_index)= ...
                norm(x(obj.local_position_coupling_indices_with_neighbor{neighbor}) - ...
                obj.x_global(ismember(obj.local_position_coupling_indices,...
                obj.local_position_coupling_indices_with_neighbor{neighbor}))')/...
                max([x(obj.local_position_coupling_indices_with_neighbor{neighbor}) ; ...
                obj.x_global(ismember(obj.local_position_coupling_indices,...
                obj.local_position_coupling_indices_with_neighbor{neighbor}))]);
            
            resid_dual_with_neighbor(rel_neig_index) =  norm(obj.rhos * ...
                (obj.x_global(ismember(obj.local_position_coupling_indices,...
                obj.local_position_coupling_indices_with_neighbor{neighbor}))- ...
                f_global_hist(ismember(obj.local_position_coupling_indices,...
                obj.local_position_coupling_indices_with_neighbor{neighbor}))))/...
                max(obj.lambdas(ismember(obj.local_position_coupling_indices,...
                obj.local_position_coupling_indices_with_neighbor{neighbor})));
            
            %set convergence allowances
            if  resid_prim_with_neighbor(rel_neig_index) < tol ...
                    && ...
                    resid_dual_with_neighbor(rel_neig_index) < tol
                convergence_allowance_to_neighbor(rel_neig_index) = 1;
            end
            
        end
        
        resid_prim = norm(x(obj.local_position_coupling_indices) - obj.x_global)/...
            max([x(obj.local_position_coupling_indices); obj.x_global]);
        resid_dual =  norm(obj.rhos * (obj.x_global-f_global_hist))/...
            max(obj.lambdas);
        
        if rem(iter,10)==0
            disp(['Iteration #',num2str(iter),': Solved cluster ', ...
                num2str(obj.cluster_names(labindex)),' (', ...
                num2str(labindex),'th cluster). total residuals', num2str(resid_prim), ', ',num2str(resid_dual),...
                '. resdls w/ neighbors: ', num2str(resid_prim_with_neighbor), ', ',num2str(resid_dual_with_neighbor),]);
        end
        
        %adjust tau
        if (1< sqrt(eta^-1*resid_prim/resid_dual)) && (sqrt(eta*resid_prim/resid_dual) < tau_max)
            tau = sqrt(eta^-1*resid_prim/resid_dual);
        elseif (tau_max^-1< sqrt(eta^-1*resid_prim/resid_dual)) && (sqrt(eta*resid_prim/resid_dual) < 1)
            tau = sqrt(eta*resid_dual/resid_prim);
        end
        
        %adjust rho
        if resid_prim > mu * resid_dual
            if obj.rhos<rho_max
                obj.rhos = obj.rhos * (tau);
            end
        elseif resid_dual > mu * resid_prim
            obj.rhos = obj.rhos / tau;
        end
        
    end
    
    
    %% check convergence
    if iter>1
        if (resid_prim < tol) && (resid_dual < tol) && (sum(convergence_allowance_from_neighbor) == numel(obj.neigh))
            converge= 1;
        else
            converge = 0;
        end
    end
    
    
    %% send flow with neighbor
    for neighbor = obj.neigh
        labSend(x(obj.local_position_coupling_indices_with_neighbor{neighbor}), find(obj.cluster_names==neighbor), maxit*labindex+iter-1);
        
        %sending the tag of the most recent data
        labSend(maxit*labindex+iter-1, find(obj.cluster_names==neighbor), maxit*100*labindex);
        
        %send convergence
        labSend(converge, find(obj.cluster_names==neighbor), maxit*200*labindex);
        
        %send convergence allowances
        labSend(convergence_allowance_to_neighbor(obj.neigh==neighbor), find(obj.cluster_names==neighbor), maxit*400*labindex);
        
    end
    
end


end
