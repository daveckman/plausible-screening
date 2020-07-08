function [SS_indicators] = ExtSTB(card_feas_region, sample_mean, sample_var, n_vec, alpha)

% Extended Screen-to-the-Best subset-selection procedure
% (flipped for minimization)

t_quantiles = tinv((1-alpha)^(1/(card_feas_region-1)), n_vec); % Quantiles for cutoffs

% Start with all solutions in subset and then eliminate
SS_indicators = ones(card_feas_region, 1);

for l = 1:card_feas_region
    for j = 1:card_feas_region
        if l ~= j
            W_lj = sqrt(t_quantiles(l)*sample_var(l)/n_vec(l) + t_quantiles(j)*sample_var(j)/n_vec(j));
            if sample_mean(l) >= sample_mean(j) + W_lj
                SS_indicators(l) = 0;
                break;
            end
        end 
    end
end