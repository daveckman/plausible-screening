function [H_QP, f_QP, A_QP, b_QP] = dCRN_formulate_QP(A, C, b, sample_mean, sample_var, n_vec)

% Take polyhedral representation of P and formulate the problem of
% minimizing the ell_CRN discrepancy as a quadratic program

% P = {(m, w): A*m + C*w <= b}
k = size(A,2);
q = size(C,2); % Number of unprojected components

% Optimization problem
% D(x_0) = min_{(m, w) in P} n*(muhat - m)'*Sigma^{-1}*(muhat - m)

H_QP = [2*nvec(1)*inv(sample_var), zeros(k, q); zeros(q, k), zeros(q, q)];
f_QP = [-2*n_vec(1)*sample_mean'/sample_var; zeros(q, 1)]; % column vector
A_QP = [A, C];
b_QP = b;

end