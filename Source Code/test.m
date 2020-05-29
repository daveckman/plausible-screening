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
crn_flag = 0;

% % TEST 1: Feasible region is {-2, -1, 0, 1, 2}^5
% feas_region = fullfact([5, 5, 5, 5, 5]) - 3;
% oracle_string = 'normacle';
% oracle_n_rngs = 1;
% exp_set = unique(unidrnd(5, [20, 5]) - 3, 'rows');
% k = size(exp_set, 1);
% n_vec = 5*ones(k, 1); % col vector
% crn_flag = 0;

%%
% SAMPLING

% Generate data and calculate summary statistics
[sample_mean, sample_var] = generate_data(oracle_string, oracle_n_rngs, exp_set, n_vec, crn_flag);

%%
% MORE SETUP
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'ell2'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'lipschitz'; % {'convex', 'lipschitz'}
prop_params = 5; % c for Lipschitz constant

%%
% SCREENING
S_indicators = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, discrep_string, fn_props, prop_params, crn_flag);
S = feas_region([1:size(feas_region,1)]'.*S_indicators,:);

%disp(S)

%S_poly_indicators = PO_screen_relax(feas_region, exp_set, sample_mean, sample_var, alpha, discrep_string, fn_props, prop_params);
%S_poly = feas_region(S_poly_indicators,:);

%%
% PLOTTING