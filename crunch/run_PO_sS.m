function run_PO_sS %(iteration, N, p)

clear;
clc;

maxNumCompThreads(4);

%add_rm_paths('add');

addAttachedFiles(gcp, {'glpk.m', 'glpkcc.mexw64'})

problem_string = 'sS_inventory';
%[oracle_string, oracle_n_rngs, feas_region, exp_set, k, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string] = init_problem(problem_string);

% (s,S) inventory problem
oracle_string = 'sS_oracle';
oracle_n_rngs = 1;
[A,B] = meshgrid(10:80,10:100);
feas_region = [A(:),B(:)];
feas_region = feas_region(feas_region(:,1)+14.5 < feas_region(:,2),:);
scrXn = [feas_region;...
    repmat(feas_region(feas_region(:,1)-feas_region(:,2)==-15,:),[5,1]);
    repmat(feas_region(feas_region(:,1)==10,:),[5,1]);
    repmat(feas_region(feas_region(:,1)==80,:),[5,1]);
    repmat(feas_region(feas_region(:,2)==10,:),[5,1]);
    repmat(feas_region(feas_region(:,2)==100,:),[5,1])];
k = 25;

% Reproduce same exp set
kmeans_rng = RandStream.create('mlfg6331_64');
opts = statset('Streams',kmeans_rng,'UseSubstreams',1);
[IDX, C] = kmeans(scrXn,25,'Options',opts);
exp_set = round(C);
exp_set(12,:) = [55,76]; % avoid singular covariance matrix my perturbing solution

n_vec = 10*ones(k, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'ell1'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'lipschitz_proj'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = 3; % gamma for Lipschitz constant
LP_solver_string = 'glpk'; %'glpk'; % {'MATLAB', 'glpk'}
clear('A', 'B', 'C', 'IDX', 'scrXn');


check_exceptions(discrep_string, fn_props, n_vec)

card_feas_region = size(feas_region, 1);

%% RUN MACROREPLICATIONS

M = 4; % Number of macroreplications

% Initialize data storage
S_indicators_d1 = zeros(card_feas_region, M);
S_indicators_d2 = zeros(card_feas_region, M);
S_indicators_dinf = zeros(card_feas_region, M);
S_poly_indicators_d1 = zeros(card_feas_region, M);
S_poly_indicators_d2 = zeros(card_feas_region, M);
S_poly_indicators_dinf = zeros(card_feas_region, M);

print_problem_header(problem_string, feas_region, exp_set, fn_props)

parfor m = 1:M
    
    
    % Sampling
    
    % Generate data using i.i.d. sampling and calculate summary statistics
    [sample_mean, sample_var] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, 'ell1');
    
    % Screening (using d1, d2, and dinf discrepancies)
    
    [S_indicators_d1(:,m), ~, S_poly_indicators_d1(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, 'ell1', fn_props, prop_params, LP_solver_string);

    [S_indicators_d2(:,m), ~, S_poly_indicators_d2(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, 'ell2', fn_props, prop_params, LP_solver_string);
    
    [S_indicators_dinf(:,m), ~, S_poly_indicators_dinf(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, 'ellinf', fn_props, prop_params, LP_solver_string);
    
    fprintf('\nRunning macrorep %d of %d.\n', m, M)
    print_screening_results('PO', 'ell1', S_indicators_d1(:,m))
    print_screening_results('PO relaxed', 'ell1', S_poly_indicators_d1(:,m))
    print_screening_results('PO', 'ell2', S_indicators_d2(:,m))
    print_screening_results('PO relaxed', 'ell2', S_poly_indicators_d2(:,m))
    print_screening_results('PO', 'ellinf', S_indicators_dinf(:,m))
    print_screening_results('PO relaxed', 'ellinf', S_poly_indicators_dinf(:,m))

end

save('sS_data2.mat')    
    
%%

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
    

% 
% %% PLOTTING SUBSETS
% 
% figure('Position', [0, 0, 600, 900])
% 
% subplot(3, 2, 1);
% plot_2d_sS_subset_prob(feas_region, exp_set, S_indicators_d1, 'PO: $d^1$')
% 
% subplot(3, 2, 2);
% plot_2d_sS_subset_prob(feas_region, exp_set, S_poly_indicators_d1, 'PO relaxed: $d^1$')
% 
% subplot(3, 2, 3);
% plot_2d_sS_subset_prob(feas_region, exp_set, S_indicators_d2, 'PO: $d^2$')
% 
% subplot(3, 2, 4);
% plot_2d_sS_subset_prob(feas_region, exp_set, S_poly_indicators_d2, 'PO relaxed: $d^2$')
% 
% subplot(3, 2, 5);
% plot_2d_sS_subset_prob(feas_region, exp_set, S_indicators_dinf, 'PO: $d^{\infty}$')
% 
% subplot(3, 2, 6);
% plot_2d_sS_subset_prob(feas_region, exp_set, S_poly_indicators_dinf, 'PO relaxed: $d^{\infty}$')
% 
% %% PLOTTING SUBSET SIZES (ECDFS)
% 
% all_S_indicators = {S_indicators_d1, S_poly_indicators_d1, S_indicators_d2, S_poly_indicators_d2, S_indicators_dinf, S_poly_indicators_dinf};
% string_names = {'PO: $d^1$', 'PO relaxed: $d^1$', 'PO: $d^2$', 'PO relaxed: $d^2$', 'PO: $d^{\infty}$', 'PO relaxed: $d^{\infty}$'};
% colors = {'b-', 'b:', 'g-', 'g:', 'm-', 'm:'};
% 
% plot_sample_size_ecdfs(all_S_indicators, string_names, colors)
% 
% %% END
% 
% add_rm_paths('remove');
