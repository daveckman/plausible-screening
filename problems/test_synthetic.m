% newsvendor problem
oracle_string = 'synthetic_grad_oracle';
oracle_n_rngs = 1;
%feas_region = [1:100]';
%exp_set = [5:10:95]';
feas_region = [0:0.01:1]';
%exp_set = [5:10:195]';
exp_set = [0.2, 0.4, 0.6, 0.8]';
k = size(exp_set, 1);

n_vec = 50*ones(k, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'ell2'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'convex'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = 0; % gamma for Lipschitz constant % = max(sell_price - cost, cost - salvage)
LP_solver_string = 'glpk'; % {'MATLAB', 'glpk'}