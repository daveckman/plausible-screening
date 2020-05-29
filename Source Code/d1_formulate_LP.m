function [f_LP, A_LP, b_LP] = d1_formulate_LP(A, C, b, sample_mean, sample_var, n_vec)

% Take polyhedral representation of P and formulate the problem of
% minimizing the ell_1 discrepancy as a linear program

% P = {(m, w): A*m + C*w <= b}
q = size(C,2); % Number of unprojected components

% Optimization problem
% D(x_0) = min_{(m, w) in P} sum_{i=1}^k sqrt{n_i}/\hat{\sigma}_i |muhat_i
% - m_i|

f_LP = [sqrt(n_vec./sample_var); sqrt(n_vec./sample_var); zeros(q, 1)]; % column vector
A_LP = [-A, A, C]; 
b_LP = b - A*sample_mean;

end

