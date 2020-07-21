function [SS_indicators] = ExtSTB(card_feas_region, sample_mean, sample_var, sample_pair_var, n_vec, alpha, discrep_string)

% Extended Screen-to-the-Best subset-selection procedure

% OR

% Screen-to-the-Best with CRN (common sample size

% (flipped for minimization)

if strcmp(discrep_string, 'CRN') == 1 % CRN
    t_quantile = tinv(1 - alpha/(card_feas_region - 1), n_vec(1) - 1); % Quantiles for cutoffs
else % i.i.d. sampling
    t_quantiles = tinv((1-alpha)^(1/(card_feas_region-1)), n_vec - 1); % Quantiles for cutoffs
end

% Start with all solutions in subset and then eliminate
SS_indicators = ones(card_feas_region, 1);

for l = 1:card_feas_region
    for j = 1:card_feas_region
        if l ~= j
            if strcmp(discrep_string, 'CRN') == 1 % CRN
                W_lj = t_quantile*sqrt(sample_pair_var(l,j)/n_vec(1));
            else % i.i.d. sampling
                W_lj = sqrt(t_quantiles(l)^2*sample_var(l)/n_vec(l) + t_quantiles(j)^2*sample_var(j)/n_vec(j));
            end
            
            if sample_mean(l) >= sample_mean(j) + W_lj
                SS_indicators(l) = 0;
                break;
            end
        end 
    end
end