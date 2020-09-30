add_rm_paths('add');

%feas_region = fullfact([3, 3, 3]) - 2;
%exp_set = feas_region([1,5,10,20],:);
%feas_region = [1; 2; 3; 4; 5];
feas_region = 3;
exp_set = [1; 2];
%exp_set = feas_region; % enumerate
k = size(exp_set, 1);
%n_vec = 30*ones(k, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'ellinf'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'convex'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = []; % gamma for Lipschitz constant
LP_solver_string = 'MATLAB';

sample_var = ones(k, 1);
n_vec = 10*ones(k, 1);
D_cutoff = calc_cutoff(k, n_vec, alpha, discrep_string);

for m = 1:10000
    if mod(m, 100) == 0
        disp(m)
    end
    sample_mean = normrnd(0, 1, [k, 1]);
    [S_indicators, D_x0s, S_poly_indicators, zs] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, discrep_string, D_cutoff, fn_props, prop_params, LP_solver_string);

    if S_indicators ~= S_poly_indicators
        print('Disagreement')
        break
    end
end