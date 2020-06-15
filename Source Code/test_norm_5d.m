% Minimize norm over discrete set in R^5
% feasible region is {-2, -1, 0, 1, 2}^5
oracle_string = 'normacle';
oracle_n_rngs = 1;
feas_region = fullfact([5, 5, 5, 5, 5]) - 3;
exp_set = unique(unidrnd(5, [10, 5]) - 3, 'rows'); % not reproducible
k = size(exp_set, 1);
n_vec = 10*ones(k, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'ell1'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'convex'; % {'convex', 'lipschitz', 'lipschitz_proj'}
prop_params = 3; % gamma for Lipschitz constant
LP_solver_string = 'glpk';