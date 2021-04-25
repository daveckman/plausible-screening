function D_cutoff  = calc_gradinftight_cutoff(k, d, n_vec, alpha)

% Calculate the uniform cutoff for the no gradients tightestdiscrepancy

n_MC_reps = 100000; % Number of Monte Carlo replications for estimation

terms = zeros(k, n_MC_reps);
for i = 1:k
    %terms(i,:) = (d+1)*(n_vec(i)-1)/(n_vec(i) - (d+1))*frnd(d + 1, n_vec(i) - (d + 1), [1, n_MC_reps]);
    %terms(i,:) = frnd(1, n_vec(i) - 1); D_cutoff = quantile(max(terms, [], 1), 1-alpha);
    terms(i,:) = trnd(n_vec(i)-1);
end
D_cutoff = (quantile(max(terms, [], 1), 1-alpha))^2;

end

