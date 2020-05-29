function S_indicators = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, discrep_string, fn_props, prop_params, crn_flag)

% Carry out plausible optima screening w.r.t. optimality.
% Use given functional properties and specified discrepancy/confidence.
% Return a column-vector of indicators on whether x_0 in S or S^c.

% Calculate cutoff
k = size(exp_set,1);
D_cutoff = calc_cutoff(k, n_vec, alpha, discrep_string);

PO_info_handle = str2func(strcat('setup_opt_',fn_props));

card_feas_region = size(feas_region, 1);
S_indicators = zeros(card_feas_region, 1);

% eventually make this a parfor
for l = 1:card_feas_region
    
    if mod(l, 10) == 0
        fprintf('Evaluating solution %d of %d', l, card_feas_region);
    end
    
    x0 = feas_region(l,:);
    
    % Setup optimization problem
    [A, C, b] = PO_info_handle(x0, exp_set, prop_params);
    
%     disp(full(A))
%     disp(full(C))
%     disp(full(b))
    
    %%% NEEED TO CONTROL WHICH DISCREPANCIES ARE ALLOWED BASED ON CRN_FLAG
    
    % Solve optimization problem
    switch discrep_string
        case 'ell1'
            % Reformulate as linear program
            [f_LP, A_LP, b_LP] = d1_formulate_LP(A, C, b, sample_mean, sample_var, n_vec);
            [~, minD_x0] = linprog(f_LP, A_LP, b_LP);
            
        case 'ell2'
            % Reformulate as quadratic program
            [H_QP, f_QP, A_QP, b_QP] = d2_formulate_QP(A, C, b, sample_mean, sample_var, n_vec);
            [~, f_val] = quadprog(H_QP, f_QP, A_QP, b_QP);
            minD_x0 = f_val + (n_vec./sample_var)'*sample_mean.^2; % Correct for constant term omitted from quadprog objective
            
        case 'ellinf'
            % Reformulate as linear program
            [f_LP, A_LP, b_LP] = dinf_formulate_LP(A, C, b, sample_mean, sample_var, n_vec);
            [~, minD_x0] = linprog(f_LP, A_LP, b_LP);
            
        case 'CRN'
            % Reformulate as quadratic program
            [H_QP, f_QP, A_QP, b_QP] = dCRN_formulate_QP(A, C, b, sample_mean, sample_var, n_vec);
            [~, f_val] = quadprog(H_QP, f_QP, A_QP, b_QP);
            minD_x0 = f_val + n_vec(1)*sample_mean'*sample_var\sample_mean; % Correct for constant term omitted from quadprog objective
            
        otherwise
            fprintf('Specify a valid discrepancy: {ell1, ell2, ellinf, CRN}.\n')
            %return 
             
    end 
    
    %disp(full(f_LP))
    %disp(full(A_LP))
    %disp(full(b_LP))
    
    disp(minD_x0)
    disp(D_cutoff)
    
    % Classify solution x0
    S_indicators(l) = (minD_x0 <= D_cutoff);
    
end
