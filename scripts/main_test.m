% Run one macroreplication of one PO algorithm on a given problem

clear;
clc;

add_rm_paths('add');

problem_string = 'tandem_budget'; % norm_k3, norm_3d, norm_5d, newsvendor, cts_newsvendor, sS_inventory, tandem_production, tandem_budget
[oracle_string, oracle_n_rngs, feas_region, exp_set, k, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string] = init_problem(problem_string);

check_exceptions(discrep_string, fn_props, n_vec)

%% CALCULATE CUTOFFS FOR PO

D_cutoff = calc_cutoff(k, n_vec, alpha, discrep_string);


%%
% SAMPLING

% Generate data and calculate summary statistics
fprintf('Generating sample data ')
if strcmp(discrep_string, 'CRN') == 1
    fprintf('with CRN... ')
else
    fprintf('with i.i.d. sampling...')
end
[sample_mean, sample_var, sample_pair_var] = generate_data(1, oracle_string, oracle_n_rngs, exp_set, n_vec, discrep_string);
fprintf('Done.\n')

%%
tic;
% SCREENING
fprintf('Screening solutions...\n')
[S_indicators, D_x0s, S_poly_indicators, zs] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, discrep_string, D_cutoff, fn_props, prop_params, LP_solver_string);
S = feas_region(S_indicators==1, :);
S_poly = feas_region(S_poly_indicators==1, :);
toc;
%%
% PRINT TO SCREEN
print_problem_header(problem_string, feas_region, exp_set, fn_props)
print_screening_results('PO', discrep_string, S_indicators)
print_screening_results('PO relax', discrep_string, S_poly_indicators)

%%
% PLOTTING
% ?

add_rm_paths('remove');