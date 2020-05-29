function [H_QP, f_QP, A_QP, b_QP] = d2_formulate_QP(A, C, b, sample_mean, sample_var, n_vec)

% Take polyhedral representation of P and formulate the problem of
% minimizing the ell_2 discrepancy as a quadratic program

% P = {(m, w): A*m + C*w <= b}
k = size(A,2);
q = size(C,2); % Number of unprojected components

% Optimization problem
% D(x_0) = min_{(m, w) in P} sum_{i=1}^k n_i/\hat{\sigma}_i^2 (muhat_i
% - m_i)^2

H_QP = [2*spdiags(n_vec./sample_var,0,k,k), zeros(k, q); zeros(q, k), zeros(q, q)];
f_QP = [-2*n_vec./sample_var.*sample_mean; zeros(q, 1)]; % column vector
A_QP = [A, C];
b_QP = b;

end