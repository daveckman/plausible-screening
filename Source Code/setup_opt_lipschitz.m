function [A, C, b] = setup_opt_lipschitz(x0, exp_set, prop_params)

% Setup optimization problem for Lipschitz performance function
% M(x_0) = proj_m(P) where P = {(m,w): A*m + C*w <= b}
% Construct A, C, and b.

% Lipschitz case
% m = (m_1, ... m_k)
% w = m_0

% m_i - m_j <= gamma*||x_i - x_j|| for all i, j = 1, ..., k with i neq j
% m_i - m_0 <= gamma*||x_i - x_0|| for all i = 1, ..., k
% -m_i + m_0 <= 0 for all i = 1, ..., k

% Lipschitz constant (gamma)
gamma = prop_params;

k = size(exp_set, 1);

num_constraints = 2*k^2;
A = spalloc(num_constraints, k, 4*k^2);
C = spalloc(num_constraints, 1, 2*k);
b = zeros(num_constraints, 1);

% Construct lhs matrices and rhs vector for Lipschitz
%%% LOOK INTO FASTER SPARSE CONSTRUCTION

% First set of constraints: 
% m_i - m_j <= gamma*||x_i - x_j|| for all i, j = 1, ..., k with i neq j

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
            b(row) = gamma*norm(exp_set(i,:) - exp_set(j,:),2);
            
            A(row + k*(k-1), i) = -1; 
            A(row + k*(k-1), j) = 1;
            b(row + k*(k-1)) = gamma*norm(exp_set(i,:) - exp_set(j,:),2);
        end
    end
end

% Second set of constraints
% m_i - m_0 <= gamma*||x_i - x_0|| for all i = 1, ..., k

A((2*k*(k-1) + 1 : 2*k*(k-1) + k),:) = speye(k);
C(2*k*(k-1) + 1 : 2*k*(k-1) + k) = -ones(k,1);
for i = 1:k
    b(2*k*(k-1) + i) = gamma*norm(exp_set(i,:) - x0);
end

% Third set of constraints
% -m_i + m_0 <= 0 for all i = 1, ..., k

A((2*k*(k-1) + k + 1 : end),:) = -speye(k);
C(2*k*(k-1) + k + 1 : end) = ones(k,1);

end

