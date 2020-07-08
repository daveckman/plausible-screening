function [A, C, b] = setup_lipschitz_proj(x0, exp_set, prop_params)

% Formulate M(x0) for Lipschitz performance function
% M(x0) = {m: A*m <= b}
% Construct A and b.

% Lipschitz case
% m = (m_1, ... m_k)

% M(x0) is described by the linear inequalities
% m_i - m_j <= gamma*min{||x_i - x_j||, ||x_i - x_0||} for all i, j = 1, ..., k with i ~= j

% Lipschitz constant (gamma)
gamma = prop_params;

k = size(exp_set, 1);

% Construct lhs matrices and rhs vector for Lipschitz

% Constraints: 
% m_i - m_j <= gamma*min{||x_i - x_j||, ||x_i - x_0||} for all i, j = 1, ..., k with i ~= j

num_constraints = k^2-k;

% Construct A
A_row = [(1:num_constraints)'; (1:num_constraints)'];
temp = repmat((1:k)', k, 1);
temp(1:k+1:end) = 0;
temp = temp(temp > 0);
A_col = [repelem((1:k)', k-1, 1); temp];
A_val = [ones(num_constraints,1); -ones(num_constraints,1)];
A = sparse(A_row, A_col, A_val, num_constraints, k);

% Construct b
temp1 = repelem(exp_set, k, 1) - repmat(exp_set, k, 1);
temp1 = sqrt(sum(temp1.^2,2));
temp2 = repelem(exp_set, k, 1) - repmat(x0, k^2, 1);
temp2 = sqrt(sum(temp2.^2,2));
temp = gamma*min(temp1, temp2);
temp(1:k+1:end) = -1;
b = temp(temp >= 0);

% Construct c
C = double.empty(num_constraints, 0);

end

