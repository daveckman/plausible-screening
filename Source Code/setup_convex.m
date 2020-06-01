function [A, C, b] = setup_convex(x0, exp_set, prop_params)

% Formulate M(x0) for convex performance function
% M(x0) = proj_m(P) where P = {(m,w): A*m + C*w <= b}
% Construct A, C, and b.

% prop_params is unused

% Convex case
% m = (m_1, ... m_k)
% w = m_0, (eta_1)^T, ..., (eta_k)^T where each eta_i in R^d

% M(x0) is described by the linear inequalities
% m_i - m_j - (eta_i)^T*(x_i - x_j) <= 0 for all i, j = 1, ..., k with i ~= j
% m_i - m_0 - (eta_i)^T*(x_i - x_0) <= 0 for all i = 1, ..., k
% -m_i + m_0 <= 0 for all i = 1, ..., k

d = size(x0,2);
k = size(exp_set, 1);

num_constraints = k^2 + k;
A = spalloc(num_constraints, k, 2*k^2);
C = spalloc(num_constraints, k*d+1, d*k^2 + 2*k);
b = zeros(num_constraints, 1);

% Construct lhs matrices and rhs vector for convex
%%% LOOK INTO FASTER SPARSE CONSTRUCTION

% First set of constraints: 
% m_i - m_j - (eta_i)^T*(x_i - x_j) <= 0 for all i, j = 1, ..., k with i ~= j

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
            C(row, 1 + (i-1)*d + 1 : 1 + i*d) = exp_set(j,:) - exp_set(i,:);
        end
    end
end

% Second set of constraints
% m_i - m_0 - (eta_i)^T*(x_i - x_0) <= 0 for all i = 1, ..., k

A((k*(k-1) + 1 : k*(k-1) + k), :) = speye(k);
C(k*(k-1) + 1 : k*(k-1) + k, 1) = -ones(k,1);
for i = 1:k
    C(k*(k-1) + i, 1 + (i-1)*d + 1 : 1 + i*d) = x0 - exp_set(i,:);
end

% Third set of constraints
% -m_i + m_0 <= 0 for all i = 1, ..., k

A((k*(k-1) + k + 1 : end), :) = -speye(k);
C(k*(k-1) + k + 1 : end, 1) = ones(k,1);

end

