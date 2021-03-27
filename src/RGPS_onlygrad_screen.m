function [Sonlygrad_d2_indicators, Sonlygrad_dinf_indicators, Sonlygrad_dinf_eff_indicators] = RGPS_onlygrad_screen(feas_region, exp_set, sample_mean, sample_mean_grad, sample_full_cov, n_vec, D_cutoff_dgrad, D_cutoff_dgradinf, D_cutoff_dgradinfeff)

% Gradient Plausible Screening with Dinfinity discrepancy
% Carry out (+ relaxed) plausible screening w.r.t. optimality and convexity.
% Return:
%   - a column-vector of indicators on whether x0 in S^poly.

k = size(exp_set, 1);
d = size(exp_set, 2);
card_feas_region = size(feas_region, 1);

Sonlygrad_d2_indicators = zeros(card_feas_region, 1);
Sonlygrad_dinf_indicators = zeros(card_feas_region, 1);
Sonlygrad_dinf_eff_indicators = zeros(card_feas_region, 1);

% RGPS with only gradients: d2 discrepancy

all_soln_conditions = zeros(k, k);
for i = 1:k
    for j = 1:k
        xdiff = exp_set(i,:) - exp_set(j,:); % x_i - x_j
        first_term = sqrt(D_cutoff_dgrad*(1/n_vec(i)*[1,-xdiff]*sample_full_cov(:,:,i)*[1;-xdiff'] + 1/n_vec(j)*sample_full_cov(1,1,j)));
        second_term = sample_mean(i) - sample_mean(j) - xdiff*sample_mean_grad(i,:)';
        all_soln_conditions(i,j) = first_term - second_term;
    end
end

if min(min(all_soln_conditions)) >= 0 

    for l = 1:card_feas_region
        %fprintf('Solution %d.\n', l)
        x0 = feas_region(l,:);
        
        % Relaxed Gradient Plausible Screening w/ only gradients
        RHS_vec = zeros(k,1);
        for i = 1:k
            xdiff = exp_set(i,:) - x0; % x_i - x0
            first_term = sqrt(D_cutoff_dgrad*(1/n_vec(i)*[0,-xdiff]*sample_full_cov(:,:,i)*[0;-xdiff']));
            second_term = xdiff*sample_mean_grad(i,:)';
            RHS_vec(i) = first_term + second_term;
        end
        
        % Classify solution x0 via relaxation
        Sonlygrad_d2_indicators(l) = (0 <= min(RHS_vec)); 

    end % end for

% otherwise all solutions are screened out
end % end if    

% RGPS with only gradients: dinfinity discrepancy

all_soln_conditions = zeros(k, k);
for i = 1:k
    for j = 1:k
        xdiff = exp_set(i,:) - exp_set(j,:); % x_i - x_j
        first_term = sqrt(D_cutoff_dgradinf*(1/n_vec(i)*[1,-xdiff]*sample_full_cov(:,:,i)*[1;-xdiff']) + sqrt(D_cutoff_dgradinf*1/n_vec(j)*sample_full_cov(1,1,j)));
        second_term = sample_mean(i) - sample_mean(j) - xdiff*sample_mean_grad(i,:)';
        all_soln_conditions(i,j) = first_term - second_term;
    end
end

if min(min(all_soln_conditions)) >= 0 

    for l = 1:card_feas_region
        %fprintf('Solution %d.\n', l)
        x0 = feas_region(l,:);
                
        % Relaxed Gradient Plausible Screening w/ only gradients
        RHS_vec = zeros(k,1);
        for i = 1:k
            xdiff = exp_set(i,:) - x0; % x_i - x0
            first_term = sqrt(D_cutoff_dgradinf*(1/n_vec(i)*[0,-xdiff]*sample_full_cov(:,:,i)*[0;-xdiff']));
            second_term = xdiff*sample_mean_grad(i,:)';
            RHS_vec(i) = first_term + second_term;
        end
        
        % Classify solution x0 via relaxation
        Sonlygrad_dinf_indicators(l) = (0 <= min(RHS_vec)); 

    end % end for

% otherwise all solutions are screened out
end % end if    

for l = 1:card_feas_region
    %fprintf('Solution %d.\n', l)
    x0 = feas_region(l,:);

    % Relaxed Gradient Plausible Screening w/ only gradients
    RHS_vec = zeros(k,1);
    for i = 1:k
        xdiff = exp_set(i,:) - x0; % x_i - x0
        first_term = sqrt(D_cutoff_dgradinfeff*(1/n_vec(i)*[0,-xdiff]*sample_full_cov(:,:,i)*[0;-xdiff']));
        second_term = xdiff*sample_mean_grad(i,:)';
        RHS_vec(i) = first_term + second_term;
    end

    % Classify solution x0 via relaxation
    Sonlygrad_dinf_eff_indicators(l) = (0 <= min(RHS_vec)); 

end % end for

end