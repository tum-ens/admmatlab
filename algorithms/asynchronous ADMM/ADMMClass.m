classdef ADMMClass
    %ADMMClass Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sub_model
        f_sub
        lin_cost_lambda
        lin_cost_rho
        diag_aux
        constraint_rows
        lambdas
        rhos
        cluster_idx
        cluster_names
        connected_lines
        neigh
        internal_plus_neighbor_var_indices
        internal_var_indices
        local_position_coupling_indices
        local_position_coupling_indices_with_neighbor
        global_int_set
        neigh_idx
        f_global
        it
        solution
        x_global
    end
    
    methods
        function obj = ADMMClass(inputArg1,inputArg2,inputArg3,inputArg4,...
                inputArg5,inputArg6, inputArg7, inputArg8, inputArg9,...
                inputArg10,inputArg11, inputArg12, inputArg13,...
                inputArg14, inputArg15, inputArg16, inputArg17)
            %ADMMClass Construct an instance of this class
            %   Detailed explanation goes here
            obj.sub_model = inputArg1;
            obj.f_sub = inputArg2;
            obj.lin_cost_lambda = inputArg3;
            obj.lin_cost_rho = inputArg4;
            obj.diag_aux = inputArg5;
            obj.constraint_rows = inputArg6;
            obj.lambdas = inputArg7;
            obj.rhos = inputArg8;
            obj.cluster_idx = inputArg9;
            obj.cluster_names = inputArg10;
            obj.neigh = inputArg11;
            obj.internal_plus_neighbor_var_indices = inputArg12;
            obj.internal_var_indices = inputArg13;
            obj.local_position_coupling_indices = inputArg14;
            obj.global_int_set = inputArg15;
            obj.neigh_idx = inputArg16;
            obj.local_position_coupling_indices_with_neighbor = inputArg17;
            obj.it = 999;
        end
       
    end
end

