function [S_poly_indicators] = PSG_screen_fast(feas_region, exp_set, sample_mean, sample_mean_grad, sample_full_cov, n_vec, D_cutoff, opttol)

% Gradient Plausible Screening with Dinfinity discrepancy
% Carry out (+ relaxed) plausible screening w.r.t. optimality and convexity.
% Return:
%   - a column-vector of indicators on whether x0 in S^poly.

k = size(exp_set, 1);
d = size(exp_set, 2);
card_feas_region = size(feas_region, 1);


S_poly_indicators = zeros(card_feas_region, 1);

for l = 1:card_feas_region
    %fprintf('Solution %d.\n', l)
    x0 = feas_region(l,:);

    % Relaxed Gradient Plausible Screening
    LHS_vec = zeros(k,1);
    for i = 1:k
        xdiff = exp_set(i,:) - x0; % x_i - x0
        first_term = -sqrt(D_cutoff*(1/n_vec(i)*[1,-xdiff]*sample_full_cov(:,:,i)*[1;-xdiff']));
        second_term = sample_mean(i) - xdiff*sample_mean_grad(i,:)';
        LHS_vec(i) = first_term + second_term;
    end
    left_max = max(LHS_vec);

    RHS_vec = zeros(k,1);
    for i = 1:k
        RHS_vec(i) = sqrt(D_cutoff*1/n_vec(i)*sample_full_cov(1,1,i)) + sample_mean(i);
    end
    right_min = min(RHS_vec) + opttol;
    
    % Second constraint to handle i=j pairs
    LHS_vec2 = zeros(k,1);
    for i = 1:k
        xdiff = exp_set(i,:) - x0; % x_i - x0
        first_term = -sqrt(D_cutoff*(1/n_vec(i)*[0,-xdiff]*sample_full_cov(:,:,i)*[0;-xdiff']));
        second_term = -xdiff*sample_mean_grad(i,:)';
        LHS_vec2(i) = first_term + second_term;
    end
    left_max2 = max(LHS_vec2);
    
    % Classify solution x0 via relaxation
    %S_poly_indicators(l) = (left_max <= right_min); 
    S_poly_indicators(l) = (left_max <= right_min && left_max2 <= opttol); 

end % end for


end