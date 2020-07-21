function run_PO_ctsnews_iid(K, N)

%clear;
clc;

fprintf('K = %d and N = %d.\n',K,N)

add_rm_paths('add');

crunch_cluster = parcluster;
%pool_obj = parpool(crunch_cluster);
maxNumCompThreads(4);
%addAttachedFiles(pool_obj, {'glpk.m', 'glpkcc.mexw64'})

% cts newsvendor problem
problem_string = 'cts_newsvendor';
oracle_string = 'cts_newsvendor_oracle';
oracle_n_rngs = 1;
%K = 20;
feas_region = [1:200]';
%exp_set = [5:10:195]';
exp_set = round((1:K)'*(200/K) - (200/(2*K)));
%k = K;
%N = 400;

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

D_cutoff_d1 = calc_cutoff(K, n_vec, alpha, 'ell1');
D_cutoff_d2 = calc_cutoff(K, n_vec, alpha, 'ell2');
D_cutoff_dinf = calc_cutoff(K, n_vec, alpha, 'ellinf');

%% RUN MACROREPLICATIONS

M = 200; % Number of macroreplications

% Initialize data storage
S_indicators_d1 = zeros(card_feas_region, M);
S_indicators_d2 = zeros(card_feas_region, M);
S_indicators_dinf = zeros(card_feas_region, M);
S_poly_indicators_d1 = zeros(card_feas_region, M);
S_poly_indicators_d2 = zeros(card_feas_region, M);
S_poly_indicators_dinf = zeros(card_feas_region, M);
SS_indicators = zeros(card_feas_region, M);

print_problem_header(problem_string, feas_region, exp_set, fn_props)

parfor (m = 1:M, crunch_cluster)

    % Generate data using i.i.d. sampling and calculate summary statistics
    [sample_mean, sample_var, ~] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, 'ell1');
    
    % Screening (using d1, d2, and dinf discrepancies)
    [S_indicators_d1(:,m), ~, S_poly_indicators_d1(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ell1', D_cutoff_d1, fn_props, prop_params, LP_solver_string);
    [S_indicators_d2(:,m), ~, S_poly_indicators_d2(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ell2', D_cutoff_d2, fn_props, prop_params, LP_solver_string);
    [S_indicators_dinf(:,m), ~, S_poly_indicators_dinf(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ellinf', D_cutoff_dinf, fn_props, prop_params, LP_solver_string);

    % Generate data using i.i.d. sampling and calculate summary statistics
    [sample_mean_SS, sample_var_SS, sample_pair_var_SS] = generate_data(m, oracle_string, oracle_n_rngs, feas_region, n_vec_SS, 'ell1');

    % Screening (using extended screen-to-the-best)
    [SS_indicators(:,m)] = ExtSTB(card_feas_region, sample_mean_SS, sample_var_SS, sample_pair_var_SS, n_vec_SS, alpha, '');
    
    fprintf('\nRunning macrorep %d of %d.\n', m, M)
    print_screening_results('PO', 'ell1', S_indicators_d1(:,m))
    print_screening_results('PO relaxed', 'ell1', S_poly_indicators_d1(:,m))
    print_screening_results('PO', 'ell2', S_indicators_d2(:,m))
    print_screening_results('PO relaxed', 'ell2', S_poly_indicators_d2(:,m))
    print_screening_results('PO', 'ellinf', S_indicators_dinf(:,m))
    print_screening_results('PO relaxed', 'ellinf', S_poly_indicators_dinf(:,m))
    print_screening_results('ESTB', '', SS_indicators(:,m))

end

save(['ctsnews_N=',num2str(N),'_K=',num2str(K),'_iid_',fn_props,'.mat'])    

add_rm_paths('remove');
