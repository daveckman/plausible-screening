% Run one macroreplication of one PO algorithm on a given problem

clear;
clc;

add_rm_paths('add');

problem_string = 'cts_newsvendor'; % norm_k3, norm_3d, norm_5d, newsvendor, cts_newsvendor, sS_inventory, tandem_production, tandem_budget
[oracle_string, oracle_n_rngs, feas_region, exp_set, k, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string] = init_problem(problem_string);

check_exceptions(discrep_string, fn_props, n_vec)
card_feas_region = size(feas_region, 1);
n_vec_SS = 20*ones(card_feas_region,1);

%% Run macroreplications

M = 200; % Number of macroreplications
SS_indicators = zeros(card_feas_region, M);

parfor m = 1:M
    
    [sample_mean_SS, sample_var_SS, sample_pair_var_SS] = generate_data(m, oracle_string, oracle_n_rngs, feas_region, n_vec_SS, 'ell1');
    [SS_indicators(:,m)] = ExtSTB(card_feas_region, sample_mean_SS, sample_var_SS, sample_pair_var_SS, n_vec_SS, alpha, '');
    print_screening_results('ESTB', '', SS_indicators(:,m))

end
    
add_rm_paths('remove');