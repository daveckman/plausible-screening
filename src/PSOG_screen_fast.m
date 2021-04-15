function [S_poly_indicators] = PSOG_screen_fast(feas_region, exp_set, sample_mean, sample_mean_grad, sample_full_cov, n_vec, D_cutoff, opttol)

% Gradient Plausible Screening with Dinfinity discrepancy
% Carry out (+ relaxed) plausible screening w.r.t. optimality and convexity.
% Return:
%   - a column-vector of indicators on whether x0 in S^poly.

k = size(exp_set, 1);
d = size(exp_set, 2);
card_feas_region = size(feas_region, 1);

S_poly_indicators = zeros(card_feas_region, 1);

% PSOG (only gradients)
for l = 1:card_feas_region
    %fprintf('Solution %d.\n', l)
    x0 = feas_region(l,:);

    LHS_vec = zeros(k, 1);
    for i = 1:k
        xdiff = exp_set(i,:) - x0; % x_i - x0
        first_term = -xdiff*sample_mean_grad(i,:)';
        second_term = -sqrt(D_cutoff*(1/n_vec(i)*[0,-xdiff]*sample_full_cov(:,:,i)*[0;-xdiff']));
        LHS_vec(i) = first_term + second_term;
    end
    
    % Classify solution x0 via relaxation
    S_poly_indicators(l) = (max(LHS_vec) <= opttol); 

end % end for
end