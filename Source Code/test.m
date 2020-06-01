% Script for testing the interface

clear;
clc;

% SETUP

% TEST 0: Toy problem
feas_region = [1; 2; 3];
oracle_string = 'normacle';
oracle_n_rngs = 1;
exp_set = [1; 3];
k = size(exp_set, 1);
n_vec = 5*ones(k, 1); % col vector

% % TEST 1: Feasible region is {-2, -1, 0, 1, 2}^5
% feas_region = fullfact([5, 5, 5, 5, 5]) - 3;
% oracle_string = 'normacle';
% oracle_n_rngs = 1;
% exp_set = unique(unidrnd(5, [20, 5]) - 3, 'rows');
% k = size(exp_set, 1);
% n_vec = 5*ones(k, 1); % col vector

%%
% MORE SETUP
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'ell2'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'lipschitz'; % {'convex', 'lipschitz'}
prop_params = 5; % c for Lipschitz constant

%%
% CHECK FOR EXCEPTIONS

accept_discreps = {'ell1', 'ell2', 'ellinf', 'CRN'};
if ~any(strcmp(discrep_string, accept_discreps))
    fprintf('ERROR: "%s" is not a valid standardized discrepancy.\n', discrep_string)
    fprintf('Please specify a valid standardized discrepancy:')
    fprintf('\t%s', accept_discreps{:})
    fprintf('.\n')
    return
end

accept_fn_props = {'lipschitz', 'convex'};
if ~any(strcmp(fn_props, accept_fn_props))
    fprintf('ERROR: "%s" is not a valid functional property.\n', fn_props)
    fprintf('Please specify a valid functional property:')
    fprintf('\t%s', accept_fn_props{:})
    fprintf('.\n')
    return
end

%%
% SAMPLING

% Generate data and calculate summary statistics
[sample_mean, sample_var] = generate_data(oracle_string, oracle_n_rngs, exp_set, n_vec, discrep_string);

%%
% SCREENING
S_indicators = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, discrep_string, fn_props, prop_params);
S = feas_region(S_indicators==1, :);

%%
% PRINT TO SCREEN

fprintf('Results of PO screening\n----------------------------------------------------\n')
fprintf('\tproblem name: \t\t\t\t\t\t%s\n', oracle_string)
fprintf('\tstandardized discrepancy: \t\t\t%s\n', discrep_string)
fprintf('\tfunctional property: \t\t\t\t%s\n', fn_props)
fprintf('\t# of solutions in feasible region: \t%d\n', size(feas_region,1))
fprintf('\t# of solutions in PO subset: \t\t%d\n', sum(S_indicators))

%%
% PLOTTING