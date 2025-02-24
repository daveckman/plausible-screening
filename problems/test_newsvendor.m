% newsvendor problem
oracle_string = 'newsvendor_oracle';
oracle_n_rngs = 1;
feas_region = [1:100]';
exp_set = [3:5:100]';
k = size(exp_set, 1);

n_vec = 25*ones(k, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'CRN'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'convex'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = 6; % gamma for Lipschitz constant % = max(sell_price - cost, cost - salvage)
LP_solver_string = 'glpk'; % {'MATLAB', 'glpk'}