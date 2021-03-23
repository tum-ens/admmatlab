function [cpx,central_soln] = runme(solver, admm_option, num_threads)

% This is the function that needs to be run by the user, for running the
% ADMM routines. 
% Three input arguments are necessary:
%
% solver (int): the ADMM approach that is employed 
%                    1: regular ADMM 
%                    2: ADMM with bin packing 
%                    3: asynchronous ADMM, 
%
% admm_option (int): the update method used in ADMM 
%                    1: Gabay (constant rho)
%                    2: Boyd  (adaptive rho)
%                    3: Song  (local rho)
%                    4: Goldstein (restart method)
%                    5: Wohlberg (relative residuals)
% 
% num_threads (int): the number of CPUs allowed by the process
%                    1: sequential run (only applicable to regular ADMM)
%                    2: parallel runs  

% The user has to set the problem matrices inside the define_input function
[A, lhs, rhs, lb, ub, f, cluster_id_from,cluster_id_to] = define_input();

select = struct();

% Manual input of ADMM parameters
select.opt.decomposition.maxiter = 2000; %default: 2000
select.opt.decomposition.mu = 10; %default: 40
select.opt.decomposition.tau = 0.1; %default: 0.05
select.opt.decomposition.Plink_init = 0e0; %default: 0e0
select.opt.decomposition.lambdas_init = 0e0; %default: 0e0
select.opt.decomposition.rhos_init = 1e2; %default: 1e2
select.opt.decomposition.convergence_tolerance = 1e-4; %default: 1e-4
select.opt.decomposition.Ccg = 10^10; %default: 10^10
select.opt.decomposition.Tf = 2; %default: 2
select.opt.decomposition.eta = 1; %default: 1
select.opt.decomposition.mu_sc = 0.999; %default: 0.999
select.opt.decomposition.iter_rho_adapt = 50; %default: 50
select.opt.decomposition.tau_max = 1; %default: 1
select.opt.decomposition.convergence_type = 1; %default: 1 (1: residuals, 2: relative cost difference)


select.opt.solver = solver;
select.opt.decomposition.admm_option = admm_option;
select.opt.num_threads=num_threads;

switch select.opt.solver
    case 1
        % Regular ADMM
        switch select.opt.decomposition.admm_option
            case 1
                if select.opt.num_threads > 1
                    disp('Solve with regular ADMM (parallel), Gabay (constant rho)')
                    [cpx,central_soln] = solve_with_ADMM_parallel_Gabay_constant_rho_release(select, A, lhs, rhs, ub, lb, f, cluster_id_from, cluster_id_to);
                else
                    disp('Solve with regular ADMM (sequential), Gabay (constant rho)')
                    [cpx,central_soln] = solve_with_ADMM_sequential_Gabay_constant_rho_release(select, A, lhs, rhs, ub, lb, f, cluster_id_from, cluster_id_to);
                end
            case 2
                if select.opt.num_threads > 1
                    disp('Solve with regular ADMM (parallel), Boyd (adaptive rho)')
                    [cpx,central_soln] = solve_with_ADMM_parallel_Boyd_adaptive_rho_release(select, A, lhs, rhs, ub, lb, f, cluster_id_from, cluster_id_to);
                else
                    disp('Solve with regular ADMM (sequential), Boyd (adaptive rho)')
                    [cpx,central_soln] = solve_with_ADMM_sequential_Boyd_adaptive_rho_release(select, A, lhs, rhs, ub, lb, f, cluster_id_from, cluster_id_to);
                end
            case 3
                if select.opt.num_threads > 1
                    disp('Solve with regular ADMM (parallel), Song (local rho)')
                    [cpx,central_soln] = solve_with_ADMM_parallel_Song_local_rho_release(select, A, lhs, rhs, ub, lb, f, cluster_id_from, cluster_id_to);
                else
                    disp('Solve with regular ADMM (sequential), Song (local rho)')
                    [cpx,central_soln] = solve_with_ADMM_sequential_Song_local_rho_release(select, A, lhs, rhs, ub, lb, f, cluster_id_from, cluster_id_to);
                end
            case 4
                if select.opt.num_threads > 1
                    disp('Solve with regular ADMM (parallel), Goldstein (restart method)')
                    [cpx,central_soln] = solve_with_ADMM_parallel_Goldstein_restart_release(select, A, lhs, rhs, ub, lb, f, cluster_id_from, cluster_id_to);
                else
                    disp('Solve with regular ADMM (sequential), Goldstein (restart method)')
                    [cpx,central_soln] = solve_with_ADMM_sequential_Goldstein_restart_release(select, A, lhs, rhs, ub, lb, f, cluster_id_from, cluster_id_to);
                end
            case 5
                if select.opt.num_threads > 1
                    disp('Solve with regular ADMM (parallel), Wohlberg (relative residuals)')
                    [cpx,central_soln] = solve_with_ADMM_parallel_Wohlberg_release(select, A, lhs, rhs, ub, lb, f, cluster_id_from, cluster_id_to);
                else
                    disp('Solve with regular ADMM (sequential), Wohlberg (relative residuals)')
                    [cpx,central_soln] = solve_with_ADMM_sequential_Wohlberg_release(select, A, lhs, rhs, ub, lb, f, cluster_id_from, cluster_id_to);
                end
        end
    case 2
        select.opt.bin_size = input('Please enter the number of CPUs allocated for each bin: (minimum: 1, maximum: available CPUs in your machine) ');
        %ADMM + Bin Packing
        switch select.opt.decomposition.admm_option
            case 1
                if select.opt.num_threads > 1
                    disp('Solve with ADMM + Bin packing (parallel), Gabay (constant rho)')
                    [cpx,central_soln] = g_solve_with_ADMM_parallel_Gabay_constant_rho_release(select, A, lhs, rhs, ub, lb, f, cluster_id_from, cluster_id_to);
                else
                    %nicht definiert f�r Bin Packing
                    error('ADMM with bin packing algorithm cannot be solved sequentially, please allow parallelism by allowing a higher number of threads (larger than 1)!')
                end
            case 2
                if select.opt.num_threads > 1
                    disp('Solve with ADMM + Bin packing (parallel), Boyd (adaptive rho)')
                    [cpx,central_soln] = g_solve_with_ADMM_parallel_Boyd_adaptive_rho_release(select, A, lhs, rhs, ub, lb, f, cluster_id_from, cluster_id_to);
                else
                    %nicht definiert f�r Bin Packing
                    error('ADMM with bin packing algorithm cannot be solved sequentially, please allow parallelism by allowing a higher number of threads (larger than 1)!')
                end
            case 3
                if select.opt.num_threads > 1
                    disp('Solve with ADMM + Bin packing (parallel), Song (local rho)')
                    [cpx,central_soln] = g_solve_with_ADMM_parallel_Song_local_rho_release(select, A, lhs, rhs, ub, lb, f, cluster_id_from, cluster_id_to);
                else
                    %nicht definiert f�r Bin Packing
                    error('ADMM with bin packing algorithm cannot be solved sequentially, please allow parallelism by allowing a higher number of threads (larger than 1)!')
                end
            case 4
                if select.opt.num_threads > 1
                    disp('Solve with ADMM + Bin packing (parallel), Goldstein (restart method)')
                    [cpx,central_soln] = g_solve_with_ADMM_parallel_Goldstein_restart_release(select, A, lhs, rhs, ub, lb, f, cluster_id_from, cluster_id_to);
                else
                    %nicht definiert f�r Bin Packing
                    error('ADMM with bin packing algorithm cannot be solved sequentially, please allow parallelism by allowing a higher number of threads (larger than 1)!')
                end
            case 5
                if select.opt.num_threads > 1
                    disp('Solve with ADMM + Bin packing (parallel), Wohlberg (relative residuals)')
                    [cpx,central_soln] = g_solve_with_ADMM_parallel_Wohlberg_release(select, A, lhs, rhs, ub, lb, f, cluster_id_from, cluster_id_to);
                else
                    %nicht definiert f�r Bin Packing
                    error('ADMM with bin packing algorithm cannot be solved sequentially, please allow parallelism by allowing a higher number of threads (larger than 1)!')
                end
            otherwise
                error('Not implemented. Please use the allowed arguments for the runme function!')
                
        end
    case 3
        %ADMM + Bin Packing + asynch
        switch select.opt.decomposition.admm_option
            case 1
                if select.opt.num_threads > 1
                    disp('Solve with asynchronous ADMM, Gabay (constant rho)')
                    [cpx,central_soln] = async_solve_with_ADMM_parallel_Boyd_constant_rho_release(select, A, lhs, rhs, ub, lb, f, cluster_id_from, cluster_id_to);
                else
                    %nicht definiert f�r Bin Packing
                    error('Asynchronous ADMM algorithm cannot be solved sequentially, please allow parallelism by allowing a higher number of threads (larger than 1)!')
                end
            otherwise
                error('Asynchronous ADMM algorithm only supports the constant rho, please set the second argument of the runme function to 1!')
                
        end
    otherwise
        error('Not implemented. Please use the allowed arguments for the runme function!')
end
