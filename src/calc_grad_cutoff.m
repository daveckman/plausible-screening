function D_cutoff  = calc_grad_cutoff(k, d, n_vec, alpha)

% Calculate the uniform cutoff for the gradient discrepancy

n_MC_reps = 100000; % Number of Monte Carlo replications for estimation

terms = zeros(k, n_MC_reps);
for i = 1:k
    % Old cutoff
    %terms(i,:) = (d+1)*(n_vec(i)-1)/(n_vec(i) - (d+1))*frnd(d + 1, n_vec(i) - (d + 1), [1, n_MC_reps]);

    % New cutoff
    terms(i,:) = 2*(n_vec(i)-1)/(n_vec(i) - 2)*frnd(2, n_vec(i) - 2, [1, n_MC_reps]);
end
D_cutoff = quantile(sum(terms, 1), 1-alpha);
        
end

