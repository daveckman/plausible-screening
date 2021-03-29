% newsvendor problem
oracle_string = 'synthetic2_grad_oracle';
oracle_n_rngs = 1;
[X, Y] = meshgrid(-1.95:.05:1.95, -1.95:.05:1.95);
n_grid = size(X,1)*size(X,2);
feas_region = [reshape(X,[n_grid, 1]), reshape(Y,[n_grid, 1])];
exp_set = [0, 0;
    -0.75, -0.75; 
    0.75, 0.75; 
    -0.75, 0.75; 
    0.75, -0.75;
    -1.5, 0;
    1.5, 0;
    0, -1.5;
    0, 1.5];
k = size(exp_set, 1);
n_vec = 20*ones(k, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'ell2'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'convex'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = 0; % gamma for Lipschitz constant % = max(sell_price - cost, cost - salvage)
LP_solver_string = 'glpk'; % {'MATLAB', 'glpk'}