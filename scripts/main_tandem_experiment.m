% Script for testing the interface

clear;
clc;

add_rm_paths('add');

problem_string = 'tandem_production';
[oracle_string, oracle_n_rngs, feas_region, exp_set, k, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string] = init_problem(problem_string);

check_exceptions(discrep_string, fn_props, n_vec)

%% CALCULATE CUTOFFS FOR PO

D_cutoff_d1 = calc_cutoff(k, n_vec, alpha, 'ell1');
D_cutoff_d2 = calc_cutoff(k, n_vec, alpha, 'ell2');
D_cutoff_dinf = calc_cutoff(k, n_vec, alpha, 'ellinf');
D_cutoff_dcrn = calc_cutoff(k, n_vec, alpha, 'CRN');

%% RUN MACROREPLICATIONS

M = 10; % Number of macroreplications
card_feas_region = size(feas_region, 1);

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
    
    % Sampling
    
    % Generate data using i.i.d. sampling and calculate summary statistics
    [sample_mean, sample_var] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, 'ell1');

    % Screening (using d1, d2, and dinf discrepancies)
    
    [S_indicators_d1(:,m), D_x0_d1, S_poly_indicators_d1(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ell1', D_cutoff_d1, fn_props, prop_params, LP_solver_string);

    [S_indicators_d2(:,m), D_x0_d2, S_poly_indicators_d2(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ell2', D_cutoff_d2, fn_props, prop_params, LP_solver_string);
    
    [S_indicators_dinf(:,m), D_x0_dinf, S_poly_indicators_dinf(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ellinf', D_cutoff_dinf, fn_props, prop_params, LP_solver_string);

    fprintf('\nRunning macrorep %d of %d.\n', m, M)
    print_screening_results('PO', 'ell1', S_indicators_d1(:,m))
    print_screening_results('PO relaxed', 'ell1', S_poly_indicators_d1(:,m))
    print_screening_results('PO', 'ell2', S_indicators_d2(:,m))
    print_screening_results('PO relaxed', 'ell2', S_poly_indicators_d2(:,m))
    print_screening_results('PO', 'ellinf', S_indicators_dinf(:,m))
    print_screening_results('PO relaxed', 'ellinf', S_poly_indicators_dinf(:,m))

end

%%
% ESTIMATE TRUE OBJECTIVE FUNCTION
M_MC = 500;
outputs = zeros(card_feas_region, M_MC);

oracle_rngs = {RandStream.create('mrg32k3a')};

for i = 1:card_feas_region
    % Extract solution x_i and sample size n_i
    x_i = feas_region(i,:);

    % Reset each stream to substream 1
    for r = 1:oracle_n_rngs
        oracle_rng = oracle_rngs{r};
        oracle_rng.Substream = 1;
    end

    % Take n_i replications at x_i
    outputs(i,:) = tandem_oracle(oracle_rngs, x_i, M_MC);
end

% Calculate summary statistics
est_true_mean = mean(outputs,2);

%% PLOTTING SUBSETS

xaxlabel = 'Mean Cycle Time of Machine 1';
yyaxlabel = 'Expected Finish Time';

figure('Position', [0, 0, 600, 900])

subplot(3, 2, 1);
plot_1d_subset_prob(feas_region(:,1), exp_set(:,1), S_indicators_d1, est_true_mean, xaxlabel, yyaxlabel, 'PO: $d^1$', alpha)

subplot(3, 2, 2);
plot_1d_subset_prob(feas_region(:,1), exp_set(:,1), S_poly_indicators_d1, est_true_mean, xaxlabel, yyaxlabel, 'PO relaxed: $d^1$', alpha)

subplot(3, 2, 3);
plot_1d_subset_prob(feas_region(:,1), exp_set(:,1), S_indicators_d2, est_true_mean, xaxlabel, yyaxlabel, 'PO: $d^2$', alpha)

subplot(3, 2, 4);
plot_1d_subset_prob(feas_region(:,1), exp_set(:,1), S_poly_indicators_d2, est_true_mean, xaxlabel, yyaxlabel, 'PO relaxed: $d^2$', alpha)

subplot(3, 2, 5);
plot_1d_subset_prob(feas_region(:,1), exp_set(:,1), S_indicators_dinf, est_true_mean, xaxlabel, yyaxlabel, 'PO: $d^{\infty}$', alpha)

subplot(3, 2, 6);
plot_1d_subset_prob(feas_region(:,1), exp_set(:,1), S_poly_indicators_dinf, est_true_mean, xaxlabel, yyaxlabel, 'PO relaxed: $d^{\infty}$', alpha)

%% PLOTTING SUBSET SIZES (ECDFS)

all_S_indicators = {S_indicators_d1, S_poly_indicators_d1, S_indicators_d2, S_poly_indicators_d2, S_indicators_dinf, S_poly_indicators_dinf};
string_names = {'PO: $d^1$', 'PO relaxed: $d^1$', 'PO: $d^2$', 'PO relaxed: $d^2$', 'PO: $d^{\infty}$', 'PO relaxed: $d^{\infty}$'};
colors = {'b-', 'b:', 'g-', 'g:', 'm-', 'm:'};

plot_sample_size_ecdfs(all_S_indicators, string_names, colors)

%% END

add_rm_paths('remove');