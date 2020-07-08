function [A, C, b] = setup_lipschitz(x0, exp_set, prop_params)

% Formulate M(x0) for Lipschitz performance function
% M(x0) = proj_m(P) where P = {(m,w): A*m + C*w <= b}
% Construct A, C, and b.

% Lipschitz case
% m = (m_1, ... m_k)
% w = m_0

% P is described by the linear inequalities
% m_i - m_j <= gamma*||x_i - x_j|| for all i, j = 1, ..., k with i ~= j
% m_i - m_0 <= gamma*||x_i - x_0|| for all i = 1, ..., k
% -m_i + m_0 <= 0 for all i = 1, ..., k

% Lipschitz constant (gamma)
gamma = prop_params;

k = size(exp_set, 1);

% Construct lhs matrices and rhs vector for Lipschitz
% For each block of constraints, construct sparse vectors describing A and C

% First set of constraints: 
% m_i - m_j <= gamma*||x_i - x_j|| for all i, j = 1, ..., k with i ~= j

first_num_constraints = k^2-k;
A_row_first = [(1:first_num_constraints)'; (1:first_num_constraints)'];
temp = repmat((1:k)', k, 1);
temp(1:k+1:end) = 0;
temp = temp(temp > 0);
A_col_first = [repelem((1:k)', k-1, 1); temp];
A_val_first = [ones(first_num_constraints,1); -ones(first_num_constraints,1)];

% C_first is all zeros

temp = repelem(exp_set, k, 1) - repmat(exp_set, k, 1);
temp = gamma*sqrt(sum(temp.^2,2));
temp(1:k+1:end) = -1;
b_first = temp(temp >= 0);

% Second set of constraints
% m_i - m_0 <= gamma*||x_i - x_0|| for all i = 1, ..., k

A_row_second = (1:k)';
A_col_second = (1:k)';
A_val_second = ones(k,1);

C_row_second = (1:k)';
C_col_second = ones(k, 1);
C_val_second = -ones(k, 1);

b_second = gamma*sqrt(sum((exp_set - repmat(x0, k, 1)).^2,2));

% Third set of constraints
% -m_i + m_0 <= 0 for all i = 1, ..., k

A_row_third = (1:k)';
A_col_third = (1:k)';
A_val_third = -ones(k,1);

C_row_third = (1:k)';
C_col_third = ones(k, 1);
C_val_third = ones(k, 1);

b_third = zeros(k, 1);

% Combine all sets of constraints
num_constraints = k^2 + k;
A = sparse([A_row_first; A_row_second + k^2 - k; A_row_third + k^2], [A_col_first; A_col_second; A_col_third], [A_val_first; A_val_second; A_val_third], num_constraints, k);
C = sparse([C_row_second + k^2 - k; C_row_third + k^2], [C_col_second; C_col_third], [C_val_second, C_val_third], num_constraints, 1);
b = [b_first; b_second; b_third];

end

