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

num_constraints = k^2-k;
A = spalloc(num_constraints, k, 4*k^2);
C = double.empty(num_constraints, 0);
b = zeros(num_constraints, 1);

% Construct lhs matrices and rhs vector for Lipschitz
%%% LOOK INTO FASTER SPARSE CONSTRUCTION

% First set of constraints: 
% m_i - m_j <= gamma*min{||x_i - x_j||, ||x_i - x_0||} for all i, j = 1, ..., k with i ~= j

for i = 1:k
    for j = 1:k
        if i ~= j            
            if j < i
                row = (i-1)*(k-1) + j;
            elseif j > i
                row = (i-1)*(k-1) + j - 1;
            end
            
            A(row, i) = 1;
            A(row, j) = -1;
            b(row) = gamma*min(norm(exp_set(i,:) - exp_set(j,:),2), norm(exp_set(i,:) - x0,2));

        end
    end
end

end

