% newsvendor problem
oracle_string = 'cts_newsvendor_oracle';
oracle_n_rngs = 1;
feas_region = [1:100]';
exp_set = [5:10:95]';
%feas_region = [1:200]';
%exp_set = [5:10:195]';
k = size(exp_set, 1);

n_vec = 20*ones(k, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'ell1'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'lipschitz_proj'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = 7; % gamma for Lipschitz constant % = max(sell_price - cost, cost - salvage)
LP_solver_string = 'glpk'; % {'MATLAB', 'glpk'}