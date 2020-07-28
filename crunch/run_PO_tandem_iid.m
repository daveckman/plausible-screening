function run_PO_tandem_iid(M, discrep_index)

%M = 1;

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
% budget = 50;
% num_machines = 5;
% 
% % # feasible solutions = (5 multichoose 50) = 316,251
% % Because of tight budget constraint (equality), reduce the dim to d = 4.
% feas_region = zeros(nchoosek(num_machines + budget - 1, budget), num_machines - 1);
% row = 1;
% for i1 = 0:budget
%     for i2 = 0:(budget - i1)
%         for i3 = 0:(budget - i1 - i2)
%             for i4 = 0:(budget - i1 - i2 - i3)
%                 feas_region(row,:) = [i1, i2, i3, i4];
%                 row = row + 1;
%             end
%         end
%     end
% end
% 
% load('tandem_budget_exp_set.mat','exp_set');
% k = size(exp_set, 1);

% Allocate budget = 25 resources across 5 machines in integral amounts
budget = 25;
num_machines = 5;

% # feasible solutions = (5 multichoose 25) = 23,751
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

% K-MEANS CONSTRUCTION OF EXP_SET
% kmeans_rng = RandStream.create('mlfg6331_64');
% opts = statset('Streams',kmeans_rng,'UseSubstreams',1);
% [~, C] = kmeans(feas_region,100,'Options',opts);
% exp_set = round(C);
% save('tandem_budget_exp_set_100.mat', 'exp_set')

load('tandem_budget_exp_set_100.mat','exp_set');
k = size(exp_set, 1);

n_vec = 100*ones(k, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_strings = {'ell1','ell2','ellinf'};
discrep_string = discrep_strings{discrep_index}; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'convex'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = []; % gamma for Lipschitz constant
LP_solver_string = 'MATLAB'; % {'MATLAB', 'glpk'}

check_exceptions(discrep_string, fn_props, n_vec)

card_feas_region = size(feas_region, 1);

%% CALCULATE CUTOFFS FOR PO

D_cutoffs = [calc_cutoff(k, n_vec, alpha, 'ell1'), calc_cutoff(k, n_vec, alpha, 'ell2'), calc_cutoff(k, n_vec, alpha, 'ellinf')];

%% RUN MACROREPLICATIONS

% Initialize data storage
S_indicators = zeros(card_feas_region, M);
S_poly_indicators = zeros(card_feas_region, M);
D_x0s = zeros(card_feas_region, M);
zs = zeros(card_feas_region, M);
PO_times = zeros(card_feas_region, M);
PO_relaxed_times = zeros(card_feas_region, M);

print_problem_header(problem_string, feas_region, exp_set, fn_props)

for m = 1:M

    % Generate data using i.i.d. sampling and calculate summary statistics
    [sample_mean, sample_var, ~] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, 'ell1');
    
    % Screening (using d1, d2, and dinf discrepancies)
    [S_indicators(:,m), D_x0s(:,m), S_poly_indicators(:,m), zs(:,m), PO_times(:,m), PO_relaxed_times(:,m)] = parPO_screen(crunch_cluster, feas_region, exp_set, sample_mean, sample_var, n_vec, discrep_string, D_cutoffs(discrep_index), fn_props, prop_params, LP_solver_string);
    
    fprintf('\nRunning macrorep %d of %d.\n', m, M)
    print_screening_results('PO', 'discrep_string', S_indicators(:,m))
    print_screening_results('PO relaxed', 'discrep_string', S_poly_indicators(:,m))

end

save(['tandem_M=',num2str(M),'_iid_',discrep_string,'_',fn_props,'.mat'])    

add_rm_paths('remove');
