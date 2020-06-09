function D_x0 = calc_min_std_discrep(discrep_string, A, C, b, sample_mean, sample_var, n_vec)

% 1. Take polyhedral representation of P.
% 2. Formulate a mathematical program for minimizing the standardized discrepancy.
% 3. Solve the mathematical program and return the minimum standardized discrepancy.

% Determine dimensions of polyhedral representation of M(x_0) = proj_m(P)
% where P = {(m, w): A*m + C*w <= b}
p = size(A,1); % Number of constraints
k = size(A,2); % Number of solutions in experimental set
q = size(C,2); % Number of unprojected components

% MATLAB's linprog() and quadprog() functions take the following arguments.

% linprog(f_LP, A_LP, b_LP, [], [], lb_LP, ub_LP) solves a linear program of the form
%   min_x f_LP'*x s.t. A_LP*x <= b_LP where lb_LP <= x <= ub_LP
%       f_LP is a column vector
%       A_LP is a matrix
%       b_LP is a column vector
%       lb_LP is a row vector
%       ub_LP is a row vector

% quadprog(h_QP, f_QP, A_QP, b_QP) solves a quadratic program of the form
%   min_x 0.5*x'*H_QP*x + f_QP'*x s.t. A_QP*x <= b_QP where
%       H_QP is a square matrix
%       f_QP is a column vector
%       A_QP is a matrix
%       b_QP is a column vector

% ** For the quadratic programs, the quadratic terms from mu_hat are omitted
% from the objective function passed to quadprog(). We add it (opt_val_offset)
% back to the optimal value returned by quadprog().**

switch discrep_string
    case 'ell1' % D_1 standardized discrepancy
            
        % Optimization problem
        % D(x_0) = min_{(m, w) in P} sum_{i=1}^k sqrt{n_i}/\hat{\sigma}_i |muhat_i - m_i|
        
        % Formulate as linear program
        f_LP = [sqrt(n_vec./sample_var); sqrt(n_vec./sample_var); zeros(q,1)];
        A_LP = [-A, A, C];
        b_LP = b - A*sample_mean;
        lb_LP = [zeros(1, 2*k), -Inf*ones(1, q)];
        ub_LP = Inf*ones(1, 2*k + q);
        
        % Solve linear prgram (suppress outputs)
        options = optimoptions('linprog','Display','none');
        [~, D_x0] = linprog(f_LP, A_LP, b_LP, [], [], lb_LP, ub_LP, options);
        
        
    case 'ell2' % D_2 standardized discrepancy
        
        % Optimization problem
        % D(x_0) = min_{(m, w) in P} sum_{i=1}^k n_i/\hat{\sigma}_i^2 (muhat_i - m_i)^2
        
        % Formulate as quadratic program  
        H_QP = [2*spdiags(n_vec./sample_var,0,k,k), zeros(k,q); zeros(q,k), zeros(q,q)];
        f_QP = [-2*n_vec./sample_var.*sample_mean; zeros(q,1)];
        A_QP = [A, C];
        b_QP = b;
        opt_val_offset = (n_vec./sample_var)'*sample_mean.^2;
        
        % Solve quadratic program (suppress outputs)
        options = optimoptions('quadprog','Display','none');
        [~, f_val] = quadprog(H_QP, f_QP, A_QP, b_QP, [], [], [], [], [], options);
        D_x0 = f_val + opt_val_offset;
        
    case 'ellinf' % D_inf standardized discrepancy

        % Optimization problem
        % D(x_0) = min_{(m, w) in P} max_{i=1}^k sqrt{n_i}/\hat{\sigma}_i |muhat_i - m_i|
        
        % Formulate as linear program
        f_LP = [zeros(k,1); zeros(q,1); 1];
        A_LP = [A, C, zeros(p, 1); -spdiags(sqrt(n_vec./sample_var),0,k,k), zeros(k,q), -ones(k,1); spdiags(sqrt(n_vec./sample_var),0,k,k), zeros(k,q), -ones(k,1)];
        b_LP = [b; -sqrt(n_vec./sample_var).*sample_mean; sqrt(n_vec./sample_var).*sample_mean];
        
        % Solve linear program (suppress outputs)
        options = optimoptions('linprog','Display','none');
        [~, D_x0] = linprog(f_LP, A_LP, b_LP, [], [], [], [], options);

    case 'CRN' % D_crn standardized discrepancy
            
        % Optimization problem
        % D(x_0) = min_{(m, w) in P} n*(muhat - m)'*Sigma^{-1}*(muhat - m)
        
        % Formulate as quadratic program
        H_QP = [2*n_vec(1)*inv(sample_var), zeros(k,q); zeros(q,k), zeros(q,q)];
        f_QP = [-2*n_vec(1)*(sample_mean'/sample_var)'; zeros(q,1)];
        A_QP = [A, C];
        b_QP = b;
        opt_val_offset = n_vec(1)*sample_mean'*(sample_var\sample_mean);
        
        % Solve quadratic program (suppress outputs)
        options = optimoptions('quadprog','Display','none');
        [~, f_val] = quadprog(H_QP, f_QP, A_QP, b_QP, [], [], [], [], [], options);
        D_x0 = f_val + opt_val_offset;
            
end % end switch

end