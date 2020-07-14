function [S_indicators, D_x0s, S_poly_indicators, zs] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, discrep_string, D_cutoff, fn_props, prop_params, LP_solver_string)

% Carry out plausible optima screening w.r.t. optimality.
% Use given functional properties and specified discrepancy/confidence.
% Return:
%   - a column-vector of indicators on whether x0 in S or S^c.
%   - a column-vector of the minimum standardized discrepancies.
%   - a column-vector of indicators on whether x0 in S^poly.
%   - a column-vector of the optimal values to the feasibility LPs.

PO_info_handle = str2func(strcat('setup_',fn_props));

card_feas_region = size(feas_region, 1);
S_indicators = zeros(card_feas_region, 1);
S_poly_indicators = zeros(card_feas_region, 1);
D_x0s = Inf*ones(card_feas_region, 1);
zs = zeros(card_feas_region, 1);

for l = 1:card_feas_region
    
    x0 = feas_region(l,:);
    
    % Setup optimization problem
    [A, C, b] = PO_info_handle(x0, exp_set, prop_params);
    
    % Calculate minimum standardized discrepancy of solution x0
    D_x0s(l) = calc_min_std_discrep(discrep_string, A, C, b, sample_mean, sample_var, n_vec, LP_solver_string);
    
    % Classify solution x0 via plausible optima approach
    S_indicators(l) = (D_x0s(l) <= D_cutoff);
    
    % Solve linear program to check feasibility
    zs(l) = check_poly_feas(discrep_string, A, C, b, sample_mean, sample_var, n_vec, D_cutoff, LP_solver_string);
    
    % Classify solution x0 via relaxation
    S_poly_indicators(l) = (zs(l) >= 0); 
        
end % end for