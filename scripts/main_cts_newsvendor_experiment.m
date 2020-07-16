% Script for testing the interface

clear;
clc;

add_rm_paths('add');

problem_string = 'cts_newsvendor';
[oracle_string, oracle_n_rngs, feas_region, exp_set, k, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string] = init_problem(problem_string);

check_exceptions(discrep_string, fn_props, n_vec)

card_feas_region = size(feas_region, 1);
per_soln_samplesize = sum(n_vec)/card_feas_region;

if floor(per_soln_samplesize) == per_soln_samplesize && per_soln_samplesize >= 2
    n_vec_SS = per_soln_samplesize*ones(card_feas_region, 1);
else
    fprintf('Total size of %d cannot be allocated equally across %d feasible solutions.\n', sum(n_vec), card_feas_region);
end

%% CALCULATE CUTOFFS FOR PO

D_cutoff_d1 = calc_cutoff(k, n_vec, alpha, 'ell1');
D_cutoff_d2 = calc_cutoff(k, n_vec, alpha, 'ell2');
D_cutoff_dinf = calc_cutoff(k, n_vec, alpha, 'ellinf');
D_cutoff_dcrn = calc_cutoff(k, n_vec, alpha, 'CRN');

%% RUN MACROREPLICATIONS

M = 200; % Number of macroreplications

% Initialize data storage
S_indicators_d1 = zeros(card_feas_region, M);
S_indicators_d2 = zeros(card_feas_region, M);
S_indicators_dinf = zeros(card_feas_region, M);
S_indicators_dcrn = zeros(card_feas_region, M);
S_poly_indicators_d1 = zeros(card_feas_region, M);
S_poly_indicators_d2 = zeros(card_feas_region, M);
S_poly_indicators_dinf = zeros(card_feas_region, M);
S_poly_indicators_dcrn = zeros(card_feas_region, M);
SS_indicators = zeros(card_feas_region, M);

print_problem_header(problem_string, feas_region, exp_set, fn_props)

parfor m = 1:M
      
%     % Sampling
%     
%     % Generate data using i.i.d. sampling and calculate summary statistics
%     [sample_mean, sample_var, ~] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, 'ell1');
%     
%     % Screening (using d1, d2, and dinf discrepancies)
%     [S_indicators_d1(:,m), D_x0_d1, S_poly_indicators_d1(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ell1', D_cutoff_d1, fn_props, prop_params, LP_solver_string);
%     [S_indicators_d2(:,m), D_x0_d2, S_poly_indicators_d2(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ell2', D_cutoff_d2, fn_props, prop_params, LP_solver_string);
%     [S_indicators_dinf(:,m), D_x0_dinf, S_poly_indicators_dinf(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ellinf', D_cutoff_dinf, fn_props, prop_params, LP_solver_string);
% 
%     fprintf('\nRunning macrorep %d of %d.\n', m, M)
%     print_screening_results('PO', 'ell1', S_indicators_d1(:,m))
%     print_screening_results('PO', 'ell2', S_indicators_d2(:,m))
%     print_screening_results('PO', 'ellinf', S_indicators_dinf(:,m))
% 
%     % 
% %     figure
% %     subplot(1,3,1)
% %     plot(D_x0_d1)
% %     subplot(1,3,2)
% %     plot(D_x0_d2)
% %     subplot(1,3,3)
% %     plot(D_x0_dinf)
% %     
    %______________________________________________________________
    
    % Generate data using i.i.d. sampling and calculate summary statistics
    [sample_mean_SS, sample_var_SS, sample_pair_var_SS] = generate_data(m, oracle_string, oracle_n_rngs, feas_region, n_vec_SS, 'ell1');

    % Screening (using extended screen-to-the-best)

    [SS_indicators(:,m)] = ExtSTB(card_feas_region, sample_mean_SS, sample_var_SS, sample_pair_var_SS, n_vec_SS, alpha, '');
    print_screening_results('ESTB', '', SS_indicators(:,m))
    
    %______________________________________________________________
% 
%     % Generate data using CRN and calculate summary statistics
%     [sample_mean, sample_var, ~] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, 'CRN');
% 
%     % Screening (using CRN discrepancy)
%     [S_indicators_dcrn(:,m), ~, S_poly_indicators_dcrn(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, 'CRN', D_cutoff_dcrn, fn_props, prop_params, LP_solver_string);
%     print_screening_results('PO', 'CRN', S_indicators_dcrn(:,m))
%      
    %______________________________________________________________

    % Generate data using CRN and calculate summary statistics
    [sample_mean_SS, sample_var_SS, sample_pair_var_SS] = generate_data(m, oracle_string, oracle_n_rngs, feas_region, n_vec_SS, 'CRN');

    % Screening (using extended screen-to-the-best)
    [SS_indicators_CRN(:,m)] = ExtSTB(card_feas_region, sample_mean_SS, sample_var_SS, sample_pair_var_SS, n_vec_SS, alpha, 'CRN');
    print_screening_results('ESTB', '', SS_indicators_CRN(:,m))
    
end

%%
% COMPUTE THE TRUE OBJECTIVE FUNCTION

% Calculate true objective function
cost = 3; % per unit order cost 
sell_price = 9; % per unit sale revenue
salvage = 1; % per unit salvage revenue
shortage = 1; % per unit shortage cost
wbl_scale = 50;
wbl_shape = 2; 

true_mean = zeros(1,length(feas_region));
neg_profit = @(D,Q) (cost*Q + shortage*(D - min(D,Q)) - sell_price*min(D,Q) - salvage*(Q - min(D,Q)));
for i = 1:length(true_mean)
    Q = feas_region(i);
    true_mean(i) = integral(@(D) neg_profit(D, Q).*wblpdf(D, wbl_scale, wbl_shape), 0, Inf);
end

%% PLOTTING SUBSETS

xaxlabel = 'Order Quantity';
yyaxlabel = 'Expected Loss';

figure('Position', [0, 0, 900, 900])

subplot(2, 2, 1);
plot_1d_subset_prob(feas_region, exp_set, SS_indicators, true_mean, xaxlabel, yyaxlabel, 'Extended Screen-to-the-Best', alpha)

subplot(2, 2, 2);
plot_1d_subset_prob(feas_region, exp_set, S_indicators_d1, true_mean, xaxlabel, yyaxlabel, 'Plausible Optima: $d^1$', alpha)

subplot(2, 2, 3);
plot_1d_subset_prob(feas_region, exp_set, S_indicators_d2, true_mean, xaxlabel, yyaxlabel, 'Plausible Optima: $d^2$', alpha)

subplot(2, 2, 4);
plot_1d_subset_prob(feas_region, exp_set, S_indicators_dinf, true_mean, xaxlabel, yyaxlabel, 'Plausible Optima: $d^{\infty}$', alpha)

%% PLOTTING SUBSET SIZES (ECDFS)

all_S_indicators = {SS_indicators, S_indicators_d1, S_indicators_d2, S_indicators_dinf};
string_names = {'ExtSTB', 'PO: $d^1$', 'PO: $d^2$', 'PO: $d^{\infty}$'};
colors = {'k:', 'b-', 'g-', 'm-'};

plot_sample_size_ecdfs(all_S_indicators, string_names, colors)

%% END

add_rm_paths('remove');
