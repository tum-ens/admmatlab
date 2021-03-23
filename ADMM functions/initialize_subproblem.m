function [sub_model, f_sub, lin_cost_lambda, lin_cost_rho, diag_aux, constraint_rows] = initialize_subproblem(select, A, lhs, rhs, ub, lb, f, cluster_idx, internal_plus_neighbor_var_indices, internal_var_indices, global_int_set)
    
    A_sub = A(full(spones(sum(abs(A(:, internal_plus_neighbor_var_indices{cluster_idx})), 2)')) == 1,:);
    constraint_rows = find(sum(abs(A_sub(:, setdiff(global_int_set,internal_var_indices{cluster_idx}))),2)==0);
    A_sub = A_sub(constraint_rows,:);
    A_sub = A_sub(:, internal_plus_neighbor_var_indices{cluster_idx});
    
    lhs_sub = lhs(full(spones(sum(abs(A(:, internal_plus_neighbor_var_indices{cluster_idx})), 2)')) == 1);
    lhs_sub = lhs_sub(constraint_rows);
    rhs_sub = rhs(full(spones(sum(abs(A(:, internal_plus_neighbor_var_indices{cluster_idx})), 2)')) == 1);
    rhs_sub = rhs_sub(constraint_rows);
    
    ub_sub = ub(internal_plus_neighbor_var_indices{cluster_idx});
    lb_sub = lb(internal_plus_neighbor_var_indices{cluster_idx});
    
    f_sub = f(internal_plus_neighbor_var_indices{cluster_idx});
    
    %den linearen Kostenterm aus lambda initializieren
    lin_cost_lambda = zeros(size(f_sub));
    
    %den linearen Kostenterm aus rho initializieren
    lin_cost_rho = zeros(size(f_sub));
    
    diag_aux = zeros(size(internal_plus_neighbor_var_indices{cluster_idx}));
    %h_sub{cluster_idx}=speye(size(internal_plus_neighbor_var_indices{cluster_idx},2));
    
    %sub_model = apply_cplex_parameters(select.opt);
    
    sub_model = struct;
    sub_model.Model.A = A_sub;
    sub_model.Model.rhs = rhs_sub;   
    sub_model.Model.lhs = lhs_sub;
    sub_model.Model.ub = ub_sub;
    sub_model.Model.lb = lb_sub;
end

