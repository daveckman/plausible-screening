function run_PO_ctsnews_crn(N, K, M)

%clear;
clc;

fprintf('N = %d and K = %d and M = %d.\n',N, K, M)

add_rm_paths('add');

crunch_cluster = parcluster;
%pool_obj = parpool(crunch_cluster);
maxNumCompThreads(4);
%addAttachedFiles(pool_obj, {'glpk.m', 'glpkcc.mexw64'})

% cts newsvendor problem
problem_string = 'cts_newsvendor';
oracle_string = 'cts_newsvendor_oracle';
oracle_n_rngs = 1;
feas_region = [1:200]';
exp_set = round((1:K)'*(200/K) - (200/(2*K)));

n_vec = (N/K)*ones(K, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'ell1'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'lipschitz_proj'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = 7; % gamma for Lipschitz constant % = max(sell_price - cost, cost - salvage)
LP_solver_string = 'MATLAB'; % {'MATLAB', 'glpk'}

check_exceptions(discrep_string, fn_props, n_vec)

card_feas_region = size(feas_region, 1);
per_soln_samplesize = sum(n_vec)/card_feas_region;

if floor(per_soln_samplesize) == per_soln_samplesize && per_soln_samplesize >= 2
    n_vec_SS = per_soln_samplesize*ones(card_feas_region, 1);
else
    fprintf('Total size of %d cannot be allocated equally across %d feasible solutions.\n', sum(n_vec), card_feas_region);
end

%% CALCULATE CUTOFFS FOR PO

D_cutoff_dcrn = calc_cutoff(K, n_vec, alpha, 'CRN');

%% RUN MACROREPLICATIONS

%M = 200; % Number of macroreplications

% Initialize data storage
S_indicators_dcrn = zeros(card_feas_region, M);
S_poly_indicators_dcrn = zeros(card_feas_region, M);
SS_indicators_CRN = zeros(card_feas_region, M);

print_problem_header(problem_string, feas_region, exp_set, fn_props)

%for m = 1:M
%parfor m = 1:M
parfor (m = 1:M, crunch_cluster)

    % Generate data using CRN and calculate summary statistics
    [sample_mean, sample_var, ~] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, 'CRN');

    % Screening (using CRN discrepancy)
    [S_indicators_dcrn(:,m), D_x0s, S_poly_indicators_dcrn(:,m), zs] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, 'CRN', D_cutoff_dcrn, fn_props, prop_params, LP_solver_string);
     
    % Generate data using CRN and calculate summary statistics
    [sample_mean_SS, sample_var_SS, sample_pair_var_SS] = generate_data(m, oracle_string, oracle_n_rngs, feas_region, n_vec_SS, 'CRN');

    % Screening (using extended screen-to-the-best)
    [SS_indicators_CRN(:,m)] = ExtSTB(card_feas_region, sample_mean_SS, sample_var_SS, sample_pair_var_SS, n_vec_SS, alpha, 'CRN');
    
    fprintf('\nRunning macrorep %d of %d.\n', m, M)
    print_screening_results('PO', 'CRN', S_indicators_dcrn(:,m))
    print_screening_results('PO (relaxed)', 'CRN', S_poly_indicators_dcrn(:,m))
    print_screening_results('ESTB', '', SS_indicators_CRN(:,m))

end

save(['ctsnews_N=',num2str(N),'_K=',num2str(K),'_M=',num2str(M),'_crn_',fn_props,'.mat'])    

add_rm_paths('remove');
