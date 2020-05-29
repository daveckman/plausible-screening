function [f_LP, A_LP, b_LP] = dinf_formulate_LP(A, C, b, sample_mean, sample_var, n_vec)

% Take polyhedral representation of P and formulate the problem of
% minimizing the ell_inf discrepancy as a linear program

% P = {(m, w): A*m + C*w <= b}
p = size(A,1);
k = size(A,2);
q = size(C,2); % Number of unprojected components

% Optimization problem
% D(x_0) = min_{(m, w) in P} max_{i=1}^k sqrt{n_i}/\hat{\sigma}_i |muhat_i
% - m_i|

f_LP = [zeros(k,1); zeros(q,1); 1]; % column vector
A_LP = [A, C, zeros(p,1); -spdiags(sqrt(n_vec./sample_var),0,k,k), zeros(k,q), -ones(k,1); spdiags(sqrt(n_vec./sample_var),0,k,k), zeros(k,q), -ones(k,1)];
b_LP = [b; -sqrt(n_vec./sample_var).*sample_mean; sqrt(n_vec./sample_var).*sample_mean];

end