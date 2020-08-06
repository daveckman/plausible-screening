function generate_tandem_data_iid(n)

%clear;
clc;

add_rm_paths('add');

crunch_cluster = parcluster;
%parpool(crunch_cluster, 'AttachedFiles', {'Par.m'});
%pool_obj = parpool(crunch_cluster);
maxNumCompThreads(8);
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

%load('tandem_budget_exp_set_100_budget50.mat','exp_set');
%k = size(exp_set, 1);

card_feas_region = size(feas_region, 1);
%n_vec = 2*ones(card_feas_region,1);

%n_vec = 100*ones(k, 1); % col vector
%alpha = 0.05; % Confidence level = 1-alpha
%discrep_strings = {'ell1','ell2','ellinf'};
%discrep_string = discrep_strings{discrep_index}; % {'ell1', 'ell2', 'ellinf', 'CRN'}
%fn_props = 'convex'; % {'convex', 'lipschitz', 'lipschitz_proj}
%prop_params = []; % gamma for Lipschitz constant
%LP_solver_string = 'MATLAB'; % {'MATLAB', 'glpk'}

%% RUN MACROREPLICATIONS

%print_problem_header(problem_string, feas_region, exp_set, fn_props)

oracle_handle = str2func(oracle_string);

% Initialize for storage
sample_mean = zeros(card_feas_region,1);
sample_var = zeros(card_feas_region,1);

tic;
m = 1;

%for i = 1:card_feas_region
parfor (i=1:card_feas_region, crunch_cluster)

    if mod(i, 100) == 0
        fprintf('Solution %d.\n', i)
    end
    
    % Extract solution x_i and sample size n_i
    x_i = feas_region(i,:);

    % Set up distinct random number streams to use
    oracle_rngs = cell(1, oracle_n_rngs);
    for r = 1:oracle_n_rngs
        oracle_rngs{r} = RandStream.create('mrg32k3a', 'NumStreams', m*oracle_n_rngs*card_feas_region, 'StreamIndices', (m - 1)*oracle_n_rngs*card_feas_region + oracle_n_rngs*(i - 1) + r);
    end

    % Take n_i replications at x_i
    outputs = oracle_handle(oracle_rngs, x_i, n);

    % Calculate summary statistics
    sample_mean(i) = mean(outputs);
    sample_var(i) = var(outputs);

end % end for
toc;

% Generate data using i.i.d. sampling and calculate summary statistics
%[sample_mean, sample_var, ~] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, 'ell1');

save(['tandem_true_mean_data_n=',num2str(n)','.mat'])    

add_rm_paths('remove');
