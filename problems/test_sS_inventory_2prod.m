% (s,S) inventory problem with 2 products
oracle_string = 'sS_2prod_oracle';
oracle_n_rngs = 2;

ss1_vec = 10:40;
QQ1_vec = 20:55;
ss2_vec = 5:35;
QQ2_vec = 5:40;

% Cross designs
card_feas_region = length(ss1_vec)*length(QQ1_vec)*length(ss2_vec)*length(QQ2_vec);
[a, b, c, d] = ndgrid(ss1_vec, QQ1_vec, ss2_vec, QQ2_vec);
feas_region = [reshape(a, card_feas_region, 1), reshape(b, card_feas_region, 1), reshape(c, card_feas_region, 1), reshape(d, card_feas_region, 1)];

ss1_vec_eval = 10:5:40;
QQ1_vec_eval = 20:5:55;
ss2_vec_eval = 5:5:35;
QQ2_vec_eval = 5:5:40;

% Cross designs
k = length(ss1_vec_eval)*length(QQ1_vec_eval)*length(ss2_vec_eval)*length(QQ2_vec_eval);
[a, b, c, d] = ndgrid(ss1_vec_eval, QQ1_vec_eval, ss2_vec_eval, QQ2_vec_eval);
exp_set = [reshape(a, k, 1), reshape(b, k, 1), reshape(c, k, 1), reshape(d, k, 1)];

n_vec = 10*ones(k, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'ell1'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'lipschitz_proj'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = 3; % gamma for Lipschitz constant
LP_solver_string = 'glpk'; % {'MATLAB', 'glpk'}
clear('card_feas_region', 'a', 'b', 'c', 'd', 'ss1_vec', 'ss1_vec_eval', 'QQ1_vec', 'QQ1_vec_eval', 'ss2_vec', 'ss2_vec_eval', 'QQ2_vec', 'QQ2_vec_eval');