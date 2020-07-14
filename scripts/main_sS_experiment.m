% Script for testing the interface

clear;
clc;

add_rm_paths('add');

problem_string = 'sS_inventory';
[oracle_string, oracle_n_rngs, feas_region, exp_set, k, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string] = init_problem(problem_string);

check_exceptions(discrep_string, fn_props, n_vec)

card_feas_region = size(feas_region, 1);

%% RUN MACROREPLICATIONS

M = 3; % Number of macroreplications

% Initialize data storage
S_indicators_d1 = zeros(card_feas_region, M);
S_indicators_d2 = zeros(card_feas_region, M);
S_indicators_dinf = zeros(card_feas_region, M);
%S_indicators_dcrn = zeros(card_feas_region, M);
S_poly_indicators_d1 = zeros(card_feas_region, M);
S_poly_indicators_d2 = zeros(card_feas_region, M);
S_poly_indicators_dinf = zeros(card_feas_region, M);
%S_poly_indicators_dcrn = zeros(card_feas_region, M);

print_problem_header(problem_string, feas_region, exp_set, fn_props)

parfor m = 1:M
    
    fprintf('\n\nRunning macrorep %d of %d.\n', m, M)
    
    % Sampling
    
    % Generate data using i.i.d. sampling and calculate summary statistics
    [sample_mean, sample_var] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, 'ell1');
    
    % Screening (using d1, d2, and dinf discrepancies)
    
    [S_indicators_d1(:,m), D_x0_d1, S_poly_indicators_d1(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, 'ell1', fn_props, prop_params, LP_solver_string);
    print_screening_results('PO', 'ell1', S_indicators_d1(:,m))
    print_screening_results('PO relaxed', 'ell1', S_poly_indicators_d1(:,m))


    [S_indicators_d2(:,m), D_x0_d2, S_poly_indicators_d2(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, 'ell2', fn_props, prop_params, LP_solver_string);
    print_screening_results('PO', 'ell2', S_indicators_d2(:,m))
    print_screening_results('PO relaxed', 'ell2', S_poly_indicators_d2(:,m))

    
    [S_indicators_dinf(:,m), D_x0_dinf, S_poly_indicators_dinf(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, 'ellinf', fn_props, prop_params, LP_solver_string);
    print_screening_results('PO', 'ellinf', S_indicators_dinf(:,m))
    print_screening_results('PO relaxed', 'ellinf', S_poly_indicators_dinf(:,m))

% 
%     figure
%     subplot(1,3,1)
%     plot(D_x0_d1)
%     subplot(1,3,2)
%     plot(D_x0_d2)
%     subplot(1,3,3)
%     plot(D_x0_dinf)
%     
    %______________________________________________________________
%     
%     % Generate data using i.i.d. sampling and calculate summary statistics
%     [sample_mean_SS, sample_var_SS] = generate_data(m, oracle_string, oracle_n_rngs, feas_region, n_vec_SS, 'ell1');
% 
%     % Screening (using extended screen-to-the-best)
% 
%     [SS_indicators(:,m)] = ExtSTB(card_feas_region, sample_mean_SS, sample_var_SS, n_vec_SS, alpha);
%     print_screening_results('ESTB', '', SS_indicators(:,m))
%     
    %______________________________________________________________


    % Generate data using CRN and calculate summary statistics
    %[sample_mean, sample_var] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, 'CRN');

    % Screening (using CRN discrepancy)
    
    %[S_indicators_dcrn(:,m), ~, S_poly_indicators_dcrn(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, 'CRN', fn_props, prop_params, LP_solver_string);
    %print_screening_results('PO', 'ellCRN', S_indicators_dinf(:,m))
     
    %______________________________________________________________

    % Generate data using CRN and calculate summary statistics
    %[sample_mean_SS, sample_var_SS] = generate_data(m, oracle_string, oracle_n_rngs, feas_region, n_vec_SS, 'ell1');

    % Screening (using extended screen-to-the-best)
    % ????
    
end

%% PLOTTING SUBSETS

figure('Position', [0, 0, 600, 900])

subplot(3, 2, 1);
plot_2d_sS_subset_prob(feas_region, exp_set, S_indicators_d1, 'PO: $d^1$')

subplot(3, 2, 2);
plot_2d_sS_subset_prob(feas_region, exp_set, S_poly_indicators_d1, 'PO relaxed: $d^1$')

subplot(3, 2, 3);
plot_2d_sS_subset_prob(feas_region, exp_set, S_indicators_d2, 'PO: $d^2$')

subplot(3, 2, 4);
plot_2d_sS_subset_prob(feas_region, exp_set, S_poly_indicators_d2, 'PO relaxed: $d^2$')

subplot(3, 2, 5);
plot_2d_sS_subset_prob(feas_region, exp_set, S_indicators_dinf, 'PO: $d^{\infty}$')

subplot(3, 2, 6);
plot_2d_sS_subset_prob(feas_region, exp_set, S_poly_indicators_dinf, 'PO relaxed: $d^{\infty}$')

%% PLOTTING SUBSET SIZES (ECDFS)

all_S_indicators = {S_indicators_d1, S_poly_indicators_d1, S_indicators_d2, S_poly_indicators_d2, S_indicators_dinf, S_poly_indicators_dinf};
string_names = {'PO: $d^1$', 'PO relaxed: $d^1$', 'PO: $d^2$', 'PO relaxed: $d^2$', 'PO: $d^{\infty}$', 'PO relaxed: $d^{\infty}$'};
colors = {'b-', 'b:', 'g-', 'g:', 'm-', 'm:'};

plot_sample_size_ecdfs(all_S_indicators, string_names, colors)

%% END

add_rm_paths('remove');
