% Run one macroreplication of one PO algorithm on a given problem

clear;
clc;

add_rm_paths('add');

problem_string = 'newsvendor'; % norm_k3, norm_3d, norm_5d, sS_inventory, tandem_production
[oracle_string, oracle_n_rngs, feas_region, exp_set, k, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string] = init_problem(problem_string);

check_exceptions(discrep_string, fn_props, n_vec)


%%
% SAMPLING

% Generate data and calculate summary statistics
fprintf('Generating sample data ')
if strcmp(discrep_string, 'CRN') == 1
    fprintf('with CRN... ')
else
    fprintf('with i.i.d. sampling...')
end
[sample_mean, sample_var] = generate_data(1, oracle_string, oracle_n_rngs, exp_set, n_vec, discrep_string);
fprintf('Done.\n')

%%
% SCREENING
fprintf('Screening solutions...\n')
[S_indicators, D_x0s, S_poly_indicators, zs] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string);
S = feas_region(S_indicators==1, :);
S_poly = feas_region(S_poly_indicators==1, :);

%%
% PRINT TO SCREEN
print_screening_results(problem_string, feas_region, exp_set, 'PO', discrep_string, fn_props, S_indicators)

%%
% PLOTTING
% ?

add_rm_paths('remove');