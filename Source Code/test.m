% Script for testing the interface

clear;
clc;

% SETUP

% % TEST 0: Toy problem
% oracle_string = 'normacle';
% oracle_n_rngs = 1;
% feas_region = [1; 2; 3];
% exp_set = [1; 3];
% k = size(exp_set, 1);
% n_vec = 5*ones(k, 1); % col vector

% % TEST 1: Feasible region is {-2, -1, 0, 1, 2}^5
% oracle_string = 'normacle';
% oracle_n_rngs = 1;
% feas_region = fullfact([5, 5, 5, 5, 5]) - 3;
% exp_set = unique(unidrnd(5, [10, 5]) - 3, 'rows');
% k = size(exp_set, 1);
% n_vec = 10*ones(k, 1); % col vector

% % TEST 2: Feasible region is {-1, 0, 1}^3
% oracle_string = 'normacle';
% oracle_n_rngs = 1;
% feas_region = fullfact([3, 3, 3]) - 2;
% exp_set = feas_region; % enumerate
% k = size(exp_set, 1);
% n_vec = 27*ones(k, 1); % col vector

% TEST 3: (s,S) inventory problem
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
K = 25;
[IDX, C] = kmeans(scrXn,25);
exp_set = round(C);
k = size(exp_set, 1);
n_vec = 10*ones(k, 1); % col vector

%%
% MORE SETUP
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'ellinf'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'lipschitz_proj'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = 3; % gamma for Lipschitz constant

%%
% CHECK FOR EXCEPTIONS

accept_discreps = {'ell1', 'ell2', 'ellinf', 'CRN'};
if ~any(strcmp(discrep_string, accept_discreps))
    fprintf('\nERROR: "%s" is not a valid standardized discrepancy.\n', discrep_string)
    fprintf('Please specify a valid standardized discrepancy:')
    fprintf('\t%s', accept_discreps{:})
    fprintf('.\n')
    return
end

accept_fn_props = {'lipschitz', 'lipschitz_proj', 'convex'};
if ~any(strcmp(fn_props, accept_fn_props))
    fprintf('\nERROR: "%s" is not a valid functional property.\n', fn_props)
    fprintf('Please specify a valid functional property:')
    fprintf('\t%s', accept_fn_props{:})
    fprintf('.\n')
    return
end

if strcmp(discrep_string, 'CRN') == 1
    % Check that all values in n_vec vector are equal
    if min(n_vec) ~= max(n_vec)
        fprintf('All sample sizes must be equal when using CRN.\n')
        return
    end
    
    % Check if n >= k so that sample_var is non-singular
    if n_vec(1) < k
        fprintf('Common sample size n = %d must be at least k = %d.\n', n_vec(1), k)
        return
    end
end

%%
% SAMPLING

% Generate data and calculate summary statistics
fprintf('Generating sample data in parallel...\n')
[sample_mean, sample_var] = generate_data(oracle_string, oracle_n_rngs, exp_set, n_vec, discrep_string);

%%
% SCREENING
fprintf('Screening solutions in parallel...\n')
[S_indicators, D_x0s, S_poly_indicators, zs] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, discrep_string, fn_props, prop_params);
S = feas_region(S_indicators==1, :);
S_poly = feas_region(S_poly_indicators==1, :);

%%
% PRINT TO SCREEN

fprintf('\nResults of PO screening\n-------------------------------------------------------\n')
fprintf('\tproblem name: \t\t\t\t\t\t\t%s\n', oracle_string)
fprintf('\tstandardized discrepancy: \t\t\t\t%s\n', discrep_string)
fprintf('\tfunctional property: \t\t\t\t\t%s\n', fn_props)
fprintf('\t# of solutions in feasible region: \t\t%d\n', size(feas_region,1))
fprintf('\t# of solutions in experimental set: \t%d\n', size(exp_set,1))
fprintf('\t# of solutions in PO subset: \t\t\t%d\n', sum(S_indicators))
fprintf('\t# of solutions in PO relaxed subset: \t%d\n', sum(S_poly_indicators))

%%
% PLOTTING