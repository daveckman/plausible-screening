function [oracle_string, oracle_n_rngs, feas_region, exp_set, k, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string] = init_problem(problem_string)
% Run a script to initialize the problem:
%   oracle_string: string for the name of the oracle function
%   oracle_n_rngs: number of rngs needed by the oracle
%   feas_region: feasible region matrix (each row corresponds to a solution)
%   exp_set: experimental set matrix (each row corresponds to a solution)
%   k: number of evaluated solutions
%   n_vec: column vector of sample sizes for each evaluated solution
%   alpha: confidence level = 1-alpha
%   discrep_string: string for discrepancy type {'ell1', 'ell2', 'ellinf', 'CRN'}
%   fn_props: string for functional property {'convex', 'lipschitz', 'lipschitz_proj}
%   prop_params: Lipschitz constant (if applicable)

%pwd
%['..\problems\test_',problem_string]
%isfile(['..\problems\test_',problem_string])
if isfile(['..\problems\test_',problem_string,'.m'])
    run(['test_',problem_string,'.m'])
else
    fprintf('\nERROR: No file with name "test_%s".\n', problem_string)
    return
end

end

