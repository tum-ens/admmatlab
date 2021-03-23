function grouping = group_subproblems_alternative(weights,bin_number)

    %Anzahl der Subprobleme (packages)
    N_p = numel(weights);
    %Anzahl der parallelen Vorgänge (bins)
    N_b = bin_number;

    grouping = cell(N_b,1);
    bin_load = zeros(N_b,1);

    weights_cp = weights;

    for i = 1:N_p
        [weight_max,idx_max] = max(weights_cp);
        [load_min,idx_min] = min(bin_load);
        bin_load(idx_min) = load_min + weight_max;
        grouping{idx_min} = [grouping{idx_min}, idx_max];
        weights_cp(idx_max) = 0;
    end
    
    for i = 1:numel(grouping)
        bin_population(i) = numel(grouping{i});
    end
    
    grouping_matrix = zeros(max(bin_population),numel(grouping));
    
    for i = 1:numel(grouping)
        for j = 1:max(bin_population(i))
            try
                grouping_matrix(j,i) = weights(grouping{i}(j));
            catch
                continue
            end
        end
    end
    bin_plot = bar(grouping_matrix','stacked')
    
end

