% Script for testing the interface

clear;
clc;

% SETUP
% Run a script to initialize:
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

% test_norm_k3;
% test_norm_3d;
% test_norm_5d;
test_sS_inventory;


%%
% CHECK FOR EXCEPTIONS

accept_discreps = {'ell1', 'ell2', 'ellinf', 'CRN'};
if ~any(strcmp(discrep_string, accept_discreps))
    fprintf('\nERROR: "%s" is not a valid standardized discrepancy.\n', discrep_string)
    fprintf('Please specify a valid standardized discrepancy:')
    fprintf('\t%s', accept_discreps{:})
    fprintf('.\n')
    return
end

accept_fn_props = {'lipschitz', 'lipschitz_proj', 'convex'};
if ~any(strcmp(fn_props, accept_fn_props))
    fprintf('\nERROR: "%s" is not a valid functional property.\n', fn_props)
    fprintf('Please specify a valid functional property:')
    fprintf('\t%s', accept_fn_props{:})
    fprintf('.\n')
    return
end

if strcmp(discrep_string, 'CRN') == 1
    % Check that all values in n_vec vector are equal
    if min(n_vec) ~= max(n_vec)
        fprintf('All sample sizes must be equal when using CRN.\n')
        return
    end
    
    % Check if n >= k so that sample_var is non-singular
    if n_vec(1) < k
        fprintf('Common sample size n = %d must be at least k = %d.\n', n_vec(1), k)
        return
    end
end

%%
% SAMPLING
m = 1; % 1 macroreplication

% Generate data and calculate summary statistics
fprintf('Generating sample data in parallel...\n')
[sample_mean, sample_var] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, discrep_string);

%%
% SCREENING
fprintf('Screening solutions in parallel...\n')
[S_indicators, D_x0s, S_poly_indicators, zs] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, discrep_string, fn_props, prop_params);
S = feas_region(S_indicators==1, :);
S_poly = feas_region(S_poly_indicators==1, :);

%%
% PRINT TO SCREEN

fprintf('\nResults of PO screening\n-------------------------------------------------------\n')
fprintf('\tproblem name: \t\t\t\t\t\t\t%s\n', oracle_string)
fprintf('\tstandardized discrepancy: \t\t\t\t%s\n', discrep_string)
fprintf('\tfunctional property: \t\t\t\t\t%s\n', fn_props)
fprintf('\t# of solutions in feasible region: \t\t%d\n', size(feas_region,1))
fprintf('\t# of solutions in experimental set: \t%d\n', size(exp_set,1))
fprintf('\t# of solutions in PO subset: \t\t\t%d\n', sum(S_indicators))
fprintf('\t# of solutions in PO relaxed subset: \t%d\n', sum(S_poly_indicators))

%%
% PLOTTING