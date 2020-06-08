% Minimize norm over discrete set in R^3
% feasible region is {-1, 0, 1}^3 (enumerated)
oracle_string = 'normacle';
oracle_n_rngs = 1;
feas_region = fullfact([3, 3, 3]) - 2;
exp_set = feas_region; % enumerate
k = size(exp_set, 1);
n_vec = 27*ones(k, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'ellinf'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'lipschitz_proj'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = 3; % gamma for Lipschitz constant