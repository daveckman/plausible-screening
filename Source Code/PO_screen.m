function [S_indicators, D_x0s] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, discrep_string, fn_props, prop_params)

% Carry out plausible optima screening w.r.t. optimality.
% Use given functional properties and specified discrepancy/confidence.
% Return a column-vector of indicators on whether x_0 in S or S^c.

% Calculate cutoff
k = size(exp_set,1);
D_cutoff = calc_cutoff(k, n_vec, alpha, discrep_string);

PO_info_handle = str2func(strcat('setup_',fn_props));

card_feas_region = size(feas_region, 1);
S_indicators = zeros(card_feas_region, 1);
D_x0s = zeros(card_feas_region, 1);

parfor_progress(card_feas_region);
parfor l = 1:card_feas_region
%for l = 1:card_feas_region
    
    x0 = feas_region(l,:);
    
    % Setup optimization problem
    [A, C, b] = PO_info_handle(x0, exp_set, prop_params);
    
    % Calculate minimum standardized discrepancy of solution x0
    D_x0s(l) = calc_min_std_discrep(discrep_string, A, C, b, sample_mean, sample_var, n_vec);
   
    % Classify solution x0
    S_indicators(l) = (D_x0s(l) <= D_cutoff);
    
    parfor_progress;
    
end % end parfor
parfor_progress(0);
