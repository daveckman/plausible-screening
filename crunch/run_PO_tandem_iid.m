function run_PO_tandem_iid(M)

%M = 1

%clear;
clc;

fprintf('M = %d.\n',M)

add_rm_paths('add');

crunch_cluster = parcluster;
%pool_obj = parpool(crunch_cluster);
maxNumCompThreads(4);
%addAttachedFiles(pool_obj, {'glpk.m', 'glpkcc.mexw64'})

% cts newsvendor problem
problem_string = 'tandem_budget';
oracle_string = 'tandem_budget_oracle';
oracle_n_rngs = 1;


% Allocate budget = 50 resources across 5 machines in integral amounts
budget = 50;
num_machines = 5;

% # feasible solutions = (5 multichoose 50) = 316,251
% Because of tight budget constraint (equality), reduce the dim to d = 4.
feas_region = zeros(nchoosek(num_machines + budget - 1, budget), num_machines - 1);
row = 1;
for i1 = 0:budget
    for i2 = 0:(budget - i1)
        for i3 = 0:(budget - i1 - i2)
            for i4 = 0:(budget - i1 - i2 - i3)
                feas_region(row,:) = [i1, i2, i3, i4];
                row = row + 1;
            end
        end
    end
end

load('tandem_budget_exp_set.mat','exp_set');
k = size(exp_set, 1);

% !!!!
%feas_region = exp_set;

n_vec = 100*ones(k, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'ellinf'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'convex'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = []; % gamma for Lipschitz constant
LP_solver_string = 'MATLAB'; % {'MATLAB', 'glpk'}

check_exceptions(discrep_string, fn_props, n_vec)

card_feas_region = size(feas_region, 1);
% per_soln_samplesize = sum(n_vec)/card_feas_region;
% 
% if floor(per_soln_samplesize) == per_soln_samplesize && per_soln_samplesize >= 2
%     n_vec_SS = per_soln_samplesize*ones(card_feas_region, 1);
% else
%     fprintf('Total size of %d cannot be allocated equally across %d feasible solutions.\n', sum(n_vec), card_feas_region);
% end

%% CALCULATE CUTOFFS FOR PO

%D_cutoff_d1 = calc_cutoff(K, n_vec, alpha, 'ell1');
%D_cutoff_d2 = calc_cutoff(K, n_vec, alpha, 'ell2');
D_cutoff_dinf = calc_cutoff(k, n_vec, alpha, 'ellinf');

%% RUN MACROREPLICATIONS

M = 1; % Number of macroreplications

% Initialize data storage
%S_indicators_d1 = zeros(card_feas_region, M);
%S_indicators_d2 = zeros(card_feas_region, M);
S_indicators_dinf = zeros(card_feas_region, M);
%S_poly_indicators_d1 = zeros(card_feas_region, M);
%S_poly_indicators_d2 = zeros(card_feas_region, M);
S_poly_indicators_dinf = zeros(card_feas_region, M);
%SS_indicators = zeros(card_feas_region, M);

print_problem_header(problem_string, feas_region, exp_set, fn_props)

for m = 1:M

    % Generate data using i.i.d. sampling and calculate summary statistics
    [sample_mean, sample_var, ~] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, 'ell1');
    
    % Screening (using d1, d2, and dinf discrepancies)
    %[S_indicators_d1(:,m), ~, S_poly_indicators_d1(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ell1', D_cutoff_d1, fn_props, prop_params, LP_solver_string);
    %[S_indicators_d2(:,m), ~, S_poly_indicators_d2(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ell2', D_cutoff_d2, fn_props, prop_params, LP_solver_string);
    [S_indicators_dinf(:), D_x0s, S_poly_indicators_dinf(:), zs] = parPO_screen(crunch_cluster, feas_region, exp_set, sample_mean, sample_var, n_vec, 'ellinf', D_cutoff_dinf, fn_props, prop_params, LP_solver_string);

    % Generate data using i.i.d. sampling and calculate summary statistics
    %[sample_mean_SS, sample_var_SS, sample_pair_var_SS] = generate_data(m, oracle_string, oracle_n_rngs, feas_region, n_vec_SS, 'ell1');

    % Screening (using extended screen-to-the-best)
    %[SS_indicators(:,m)] = ExtSTB(card_feas_region, sample_mean_SS, sample_var_SS, sample_pair_var_SS, n_vec_SS, alpha, '');
    
    fprintf('\nRunning macrorep %d of %d.\n', m, M)
    %print_screening_results('PO', 'ell1', S_indicators_d1(:,m))
    %print_screening_results('PO relaxed', 'ell1', S_poly_indicators_d1(:,m))
    %print_screening_results('PO', 'ell2', S_indicators_d2(:,m))
    %print_screening_results('PO relaxed', 'ell2', S_poly_indicators_d2(:,m))
    print_screening_results('PO', 'ellinf', S_indicators_dinf(:,m))
    print_screening_results('PO relaxed', 'ellinf', S_poly_indicators_dinf(:,m))
    %print_screening_results('ESTB', '', SS_indicators(:,m))

end

save(['tandem_M=',num2str(M),'_iid_',fn_props,'.mat'])    

add_rm_paths('remove');
