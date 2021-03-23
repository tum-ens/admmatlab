function grouping = group_subproblems(weights,bin_number)
%clear
addpath(fullfile('C:\gurobi811', 'win64', 'examples', 'matlab'));

%Anzahl der Subprobleme (packages)
N_p = numel(weights);
package_series= [1:N_p];

%Anzahl der parallelen Vorgänge (bins)
N_b = bin_number;
bin_series = [1:N_b];

weights_variable = weights';
weights_original = weights_variable;

S = optimvar('size');
record_soln= [];

x = optimvar('x', N_b, N_p,'Type','integer','LowerBound',0,'UpperBound',1);

prob = cell(N_b,1);
soln = cell(N_b,1);
prob{1} = optimproblem('ObjectiveSense','minimize');

prob{1}.Constraints.all_problems_packed = [sum(x) == ones(1,N_p)];
prob{1}.Constraints.bin_capacity_respected = [x * weights_variable <= S * ones(N_b,1)];

prob{1}.Objective = S;
options = optimoptions('intlinprog','RelativeGapTolerance',1e-3);
soln{1}=solve(prob{1},'Options',options);

record_soln = 999*ones(N_b,N_p);
%save solution
[~, max_bin_local_index] = max(soln{1}.x*weights_variable);
remove_package_local_index = find(soln{1}.x(max_bin_local_index,:));
bins_left = 1:N_b;

record_soln(bins_left,remove_package_local_index) = soln{1}.x(:,remove_package_local_index);
remove_package_global_index = remove_package_local_index;
max_bin_global_index = max_bin_local_index;
record_soln
bins_left(max_bin_global_index) = []; 
N_b_org=N_b;
for iter=1:N_b_org-1
    package_mapping=[];
    %identify maximum bin
    %[~, max_bin_local_index] = max(soln{iter}.x*w_p);
    %identify the packages in the maximum bin
    %remove_package_local_index = find(soln{iter}.x(max_bin_local_index,:));

    %identify package mapping for next iteration
    package_series(remove_package_local_index) = [];
    package_mapping(1:size(soln{iter}.x,2)-numel(remove_package_local_index)) = package_series;

    %identify bin mapping for next iteration
    bin_series(max_bin_local_index) = [];
    bin_mapping(1:size(soln{iter}.x,1)-1) = bin_series;
    
    %remove the maximum bin from the last iteration
    %solve again
    prob{iter+1} = optimproblem('ObjectiveSense','minimize');
    
    N_b = N_b - 1;
    N_p = N_p - numel(remove_package_global_index);
    
    x = optimvar('x', N_b, N_p,'Type','integer','LowerBound',0,'UpperBound',1);
    weights_variable(remove_package_local_index)=[];
    prob{iter+1}.Constraints.all_problems_packed = [sum(x,1) == ones(1,N_p)];
    prob{iter+1}.Constraints.bin_capacity_respected = [x * weights_variable <= S * ones(N_b,1)];   
    prob{iter+1}.Objective = S;

    soln{iter+1}=solve(prob{iter+1},'Options',options);
    
    [~, max_bin_local_index] = max(soln{iter+1}.x*weights_variable);
    remove_package_local_index = find(soln{iter+1}.x(max_bin_local_index,:));
    
    max_bin_global_index = bin_mapping(max_bin_local_index);
    
    remove_package_global_index = package_mapping(remove_package_local_index);
    
    record_soln(bins_left,remove_package_global_index) = soln{iter+1}.x(:,remove_package_local_index);
    bins_left = bins_left(bins_left~=max_bin_global_index); 

end
record_soln(record_soln==999) = 0;
solution = record_soln*weights_original;

grouping=cell(N_b_org,1);
for iter=1:N_b_org
    grouping{iter} = find(record_soln(iter,:));
end

