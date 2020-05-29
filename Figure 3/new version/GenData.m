Xr = repmat(X,[num_sim ,1]);
Xu = unique(X, 'rows');

muhat_by_macro = zeros(M, K);
sigma2hat_by_macro = zeros(M, K);

for m = 1:M
    raw_data = sSsimu(Xr);
    
    for k = 1:K
        Y = raw_data(k:K:num_sim*K);
        muhat_by_macro(m,k) = mean(Y);
        sigma2hat_by_macro(m,k) = var(Y);
    end
end