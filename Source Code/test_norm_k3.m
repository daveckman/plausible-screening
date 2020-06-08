% Toy problem
oracle_string = 'normacle';
oracle_n_rngs = 1;
feas_region = [1; 2; 3];
exp_set = [1; 3];
k = size(exp_set, 1);
n_vec = 5*ones(k, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'ellinf'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'lipschitz_proj'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = 3; % gamma for Lipschitz constant