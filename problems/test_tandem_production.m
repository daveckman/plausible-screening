% tandem production line problem
oracle_string = 'tandem_oracle';
oracle_n_rngs = 1;
%feas_region = [0.5, 0.5; 0.4, 0.6];
feas_region = [(1:99)', flip(1:99)']./100;
exp_set = feas_region(4:9:99,:);
k = size(exp_set, 1);

n_vec = 30*ones(k, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'CRN'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'convex'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = []; % gamma for Lipschitz constant
LP_solver_string = 'glpk'; % {'MATLAB', 'glpk'}