function [A, C, b] = setup_convex_grad(x0, exp_set)

% Formulate Z(x0) for convex performance function
% Z(x0) = proj_z(P) where Z = {(z,w): A*z + C*w <= b}
% Construct A, C, and b.

% prop_params is unused

% Convex case
% z = (m_1, (g_11, ..., g_1d), ..., m_k, (g_k1, ..., g_kd))
% w = m_0

% P is described by the linear inequalities
% m_i - m_j - (x_i - x_j)^T*g_i <= 0 for all i, j = 1, ..., k with i ~= j
% m_i - m_0 - (x_i - x_0)^T*g_i <= 0 for all i = 1, ..., k
% -m_i + m_0 <= 0 for all i = 1, ..., k

d = size(x0,2);
k = size(exp_set, 1);

% Construct lhs matrices and rhs vector for convex
% For each block of constraints, construct sparse vectors describing A and C

% First set of constraints: 
% m_i - m_j - (eta_i)^T*(x_i - x_j) <= 0 for all i, j = 1, ..., k with i ~= j

first_num_constraints = k^2-k;
A_row_first = [(1:first_num_constraints)'; (1:first_num_constraints)'];
temp = repmat((1:k)', k, 1);
temp(1:k+1:end) = 0;
temp = temp(temp > 0);
A_col_first = [repelem((1:k)', k-1, 1); temp];
A_val_first = [ones(first_num_constraints,1); -ones(first_num_constraints,1)];

C_row_first = repelem((1:(k^2-k))', d, 1);
temp1 = repmat(1:d,(k-1)*k, 1);
temp2 = repmat(repelem((0:k-1)',k-1,1), 1, d);
temp3 = d*temp2 + temp1;
C_col_first = 1 + reshape(temp3', d*(k-1)*k, 1);

temp = repmat(exp_set, k, 1) - repelem(exp_set, k, 1);
temp(1:k+1:end,:) = [];
C_val_first = reshape(temp', d*(k-1)*k, 1);

% Second set of constraints
% m_i - m_0 - (eta_i)^T*(x_i - x_0) <= 0 for all i = 1, ..., k

A_row_second = (1:k)';
A_col_second = (1:k)';
A_val_second = ones(k,1);

% m0 terms
C_row_second_m0 = (1:k)';
C_col_second_m0 = ones(k,1);
C_val_second_m0 = -ones(k,1);

% eta terms
C_row_second_etas = repelem((1:k)', d, 1);
C_col_second_etas = (1:k*d)' + 1;
C_val_second_etas = repmat(x0', k, 1) - reshape(exp_set', k*d, 1);

% Combine C vectors
C_row_second = [C_row_second_m0; C_row_second_etas];
C_col_second = [C_col_second_m0; C_col_second_etas];
C_val_second = [C_val_second_m0; C_val_second_etas];

% Third set of constraints
% -m_i + m_0 <= 0 for all i = 1, ..., k

A_row_third = (1:k)';
A_col_third = (1:k)';
A_val_third = -ones(k,1);

C_row_third = (1:k)';
C_col_third = ones(k, 1);
C_val_third = ones(k, 1);

% Combine all sets of constraints
num_constraints = k^2 + k;
A = sparse([A_row_first; A_row_second + k^2 - k; A_row_third + k^2], [A_col_first; A_col_second; A_col_third], [A_val_first; A_val_second; A_val_third], num_constraints, k);
C = sparse([C_row_first; C_row_second + k^2 - k; C_row_third + k^2], [C_col_first; C_col_second; C_col_third], [C_val_first; C_val_second; C_val_third], num_constraints, 1+k*d);
b = zeros(num_constraints, 1);

% FOR GRADIENTS:
% move dummy variables over to A
A = [A, C(:,1:k*d)];
% reorder columns of A
concat = [1:k; reshape(k+(1:(k*d)), d, k)];
concat = reshape(concat, k*(d+1), 1);
A = A(:,concat);
C = C(:,k*d+1);

end

