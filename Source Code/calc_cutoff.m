function D_cutoff  = calc_cutoff(k, n_vec, alpha, discrep_string)

% Calculate the uniform cutoff for a given standardized discrepancy

n_MC_reps = 10000; % Number of Monte Carlo replications for estimation

switch discrep_string
    case 'ell1'
        terms = zeros(k, n_MC_reps);
        for i = 1:k
            terms(i,:) = trnd(n_vec(i)-1, [1, n_MC_reps]);
        end
        D_cutoff = quantile(sum(abs(terms), 1), 1-alpha);
        
    case 'ell2'
        terms = zeros(k, n_MC_reps);
        for i = 1:k
            terms(i,:) = frnd(1, n_vec(i) - 1, [1, n_MC_reps]);
        end
        D_cutoff = quantile(sum(terms, 1), 1-alpha);
        
    case 'ellinf'
        terms = zeros(k, n_MC_reps);
        for i = 1:k
            terms(i,:) = trnd(n_vec(i)-1, [1, n_MC_reps]);
        end
        D_cutoff = quantile(max(abs(terms), [], 1), 1-alpha);
        
    case 'CRN'
        D_cutoff = k*(n_vec(1) - 1)/(n_vec(1) - k) * finv(1-alpha, k, n_vec(1)-k);
        
    otherwise
        fprintf('Specify a valid discrepancy: {ell1, ell2, ellinf, CRN}.\n')
        return
end

end
