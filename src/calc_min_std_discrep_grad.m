function D_x0 = calc_min_std_discrep_grad(A, C, b, sample_mean, sample_mean_grad, sample_full_cov, n_vec)

% 1. Take polyhedral representation of P.
% 2. Formulate a mathematical program for minimizing the standardized discrepancy.
% 3. Solve the mathematical program and return the minimum standardized discrepancy.

% Determine dimensions of polyhedral representation of M(x_0) = proj_m(P)
% where P = {(m, w): A*m + C*w <= b}
p = size(A,1); % Number of constraints
q = size(C,2); % Number of unprojected components
k = size(sample_mean_grad,1); % Number of solutions in experimental set
d = size(sample_mean_grad,2); % Dimension of solution space

% Concatenate sample means and gradients
concat = [sample_mean'; sample_mean_grad'];
zetahat = reshape(concat, k*(d+1), 1);

% quadprog(h_QP, f_QP, A_QP, b_QP) solves a quadratic program of the form
%   min_x 0.5*x'*H_QP*x + f_QP'*x s.t. A_QP*x <= b_QP where
%       H_QP is a square matrix
%       f_QP is a column vector
%       A_QP is a matrix
%       b_QP is a column vector

% Gradient Quadratic Program

% Optimization problem
% D(x_0) = min_{(z, w) in P} sum_i n_i*(zetahat_i - z_i)'*Psi_i^{-1}*(zetaha_i - z_i)

% Add correct to diagonal before taking inverse: epsilon = 1
%epsilon = 1;

H_QP = sparse(k*(d+1)+q, k*(d+1)+q); % with padding
for i = 1:k
    %H_QP((d+1)*(i-1)+1 : (d+1)*i, (d+1)*(i-1)+1 : (d+1)*i) = 0.5*n_vec(i)*inv(sample_full_cov(:,:,i) + epsilon*eye(d+1));
    H_QP((d+1)*(i-1)+1 : (d+1)*i, (d+1)*(i-1)+1 : (d+1)*i) = 0.5*n_vec(i)*inv(sample_full_cov(:,:,i)); %+ epsilon*eye(d+1));
end
% Add another correction
H_QP = H_QP + 0.01*speye(k*(d+1)+q);

%%

f_QP = sparse(k*(d+1)+q, 1);
A_QP = [-A, C];
b_QP = b - A*zetahat;

options = optimoptions('quadprog','Display','none','MaxIter',500); %,'OptimalityTolerance',10^(-3)); % Default tolerance = 1e-8
[~, D_x0, exitflag] = quadprog(H_QP, f_QP, A_QP, b_QP, [], [], [], [], [], options);
%display(exitflag)

% % Check for solutions with sample variances of zero.
% if strcmp(discrep_string, 'CRN') == 1
%     zero_var_solns = (diag(sample_var) == 0)'; % vector of indicators for zero-variance solutions
% else
%     zero_var_solns = (sample_var == 0)'; % vector of indicators for zero-variance solutions
% end
% k_bar = k - sum(zero_var_solns); % number of solutions with non-zero variance
% if k_bar ~= k % If any, set m_i = mu_hat_i and move from A (LHS) to b (RHS).
%     b = b - A(:,zero_var_solns)*sample_mean(zero_var_solns);
%     A(:,zero_var_solns) = [];
% end

% % MATLAB's linprog() and quadprog() functions take the following arguments.
% if nargin > 8
%     init_sol = varargin{1};
% else
%     init_sol = [];
% end

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

% 
% % Optimization problem
% % D(x_0) = min_{(z, w) in P} sum_i n_i*(zetahat_i - z_i)'*Psi_i^{-1}*(zetaha_i - z_i)
% 
% % Formulate as quadratic program
% %sample_var = sample_var + 0.01*speye(k); % Crude regularization
% H_QP = [2*n_vec(1)*fullCovinv, zeros(k,q); zeros(q,k), zeros(q,q)];
% f_QP = [-2*n_vec(1)*(sample_mean(~zero_var_solns)'/sample_var(~zero_var_solns,~zero_var_solns))'; zeros(q,1)];
% A_QP = [A, C];
% b_QP = b;
% opt_val_offset = n_vec(1)*sample_mean(~zero_var_solns)'*(sample_var(~zero_var_solns,~zero_var_solns)\sample_mean(~zero_var_solns));
% 
% % Solve quadratic program (suppress outputs) %'Display','none',
% options = optimoptions('quadprog','Display','none','MaxIter',500); %,'OptimalityTolerance',10^(-3)); % Default tolerance = 1e-8
% [~, f_val, exitflag] = quadprog(H_QP, f_QP, A_QP, b_QP, [], [], [], [], [], options);
% D_x0 = f_val + opt_val_offset;
% %         if exitflag ~= 1
% %             fprintf('\nExit flag = %d. Outputting %f.',exitflag, D_x0)
% %         end
% 
% if exitflag == -2 % infeasible QP
%     D_x0 = Inf;
% end
 

% switch discrep_string
%     case 'ell1' % D_1 standardized discrepancy
%         
%         
%         % Optimization problem
%         % D(x_0) = min_{(m, w) in P} sum_{i=1}^k sqrt{n_i}/\hat{\sigma}_i |muhat_i - m_i|
%         
%         % Formulate as linear program
%         f_LP = [sqrt(n_vec(~zero_var_solns)./sample_var(~zero_var_solns)); sqrt(n_vec(~zero_var_solns)./sample_var(~zero_var_solns)); zeros(q,1)];
%         A_LP = [-A, A, C];
%         b_LP = b - A*sample_mean(~zero_var_solns);
%         lb_LP = [zeros(1, 2*k_bar), -Inf*ones(1, q)];
%         ub_LP = Inf*ones(1, 2*k_bar + q);
%         
%         % Solve linear prgram (suppress outputs)
%         switch LP_solver_string
%             case 'MATLAB'        
%                 options = optimoptions('linprog','Display','none'); %,'OptimalityTolerance',10^(-3)); % Default tolerance = 1e-8
%                 [~, D_x0, exitflag] = linprog(f_LP, A_LP, b_LP, [], [], lb_LP, ub_LP, options);
%                 if exitflag == -2 % infeasible LP
%                     D_x0 = Inf;
%                 end
%             case 'glpk'        
%                 [~, D_x0, status] = glpkcc(f_LP, A_LP, b_LP, lb_LP, ub_LP, repmat('U',size(A_LP,1),1), repmat('C',size(A_LP,2),1), 1, struct('savefilename','SimpleLP'));
%                 if status == 110 % infeasible LP
%                     D_x0 = Inf;
%                 end
%         end
%         
%     case 'ell2' % D_2 standardized discrepancy
%         
%         % Optimization problem
%         % D(x_0) = min_{(m, w) in P} sum_{i=1}^k n_i/\hat{\sigma}_i^2 (muhat_i - m_i)^2
%         
%         % Formulate as quadratic program  
%         H_QP = [2*spdiags(n_vec(~zero_var_solns)./sample_var(~zero_var_solns),0,k_bar,k_bar), zeros(k_bar,q); zeros(q,k_bar), zeros(q,q)];
%         H_QP_reg = H_QP + 0.00001*speye(k_bar + q);
%         f_QP = [-2*n_vec(~zero_var_solns)./sample_var(~zero_var_solns).*sample_mean(~zero_var_solns); zeros(q,1)];
%         A_QP = [A, C];
%         b_QP = b;
%         opt_val_offset = (n_vec(~zero_var_solns)./sample_var(~zero_var_solns))'*sample_mean(~zero_var_solns).^2;
%         
%         % Solve quadratic program (suppress outputs)
%         options = optimoptions('quadprog','Display','none'); %,'OptimalityTolerance',10^(-3)); % Default tolerance = 1e-8
%         %init_sol = [((n_vec./sample_var)'*sample_mean)/(n_vec'*sample_var)*ones(k,1); zeros(q, 1)];
%         [~, f_val, exitflag] = quadprog(H_QP_reg, f_QP, A_QP, b_QP, [], [], [], [], init_sol, options);
%         D_x0 = f_val + opt_val_offset;
%         if exitflag == -2 % infeasible QP
%             D_x0 = Inf;
%         end
%         
%         % GLPK options (slower and less reliable)       
%         %[~, f_val] = qpng(H_QP_reg, f_QP, A_QP, b_QP);
%         %init_sol = [((n_vec./sample_var)'*sample_mean)/(n_vec'*sample_var)*ones(k,1); zeros(q, 1)];
%         %[x_opt, ~, ~, ~] = qpsolng(H_QP_reg, f_QP, [], [], A_QP, b_QP, init_sol);
%         %f_val = 0.5 * x_opt' * H_QP * x_opt + f_QP' * x_opt; 
%         %D_x0 = f_val + opt_val_offset;
%         
%     case 'ellinf' % D_inf standardized discrepancy
% 
%         % Optimization problem
%         % D(x_0) = min_{(m, w) in P} max_{i=1}^k sqrt{n_i}/\hat{\sigma}_i |muhat_i - m_i|
%         
%         % Formulate as linear program
%         f_LP = [zeros(k_bar,1); zeros(q,1); 1];
%         A_LP = [A, C, zeros(p, 1); -spdiags(sqrt(n_vec(~zero_var_solns)./sample_var(~zero_var_solns)),0,k_bar,k_bar), zeros(k_bar,q), -ones(k_bar,1); spdiags(sqrt(n_vec(~zero_var_solns)./sample_var(~zero_var_solns)),0,k_bar,k_bar), zeros(k_bar,q), -ones(k_bar,1)];
%         b_LP = [b; -sqrt(n_vec(~zero_var_solns)./sample_var(~zero_var_solns)).*sample_mean(~zero_var_solns); sqrt(n_vec(~zero_var_solns)./sample_var(~zero_var_solns)).*sample_mean(~zero_var_solns)];
%         
%         % Solve linear program (suppress outputs)
%         switch LP_solver_string
%             case 'MATLAB'        
%                 options = optimoptions('linprog','Display','none'); %,'OptimalityTolerance',10^(-3)); % Default tolerance = 1e-8
%                 [~, D_x0, exitflag] = linprog(f_LP, A_LP, b_LP, [], [], [], [], options);
%                 if exitflag == -2 % infeasible LP
%                     D_x0 = Inf;
%                 end
%             case 'glpk'
%                 [~, D_x0, status] = glpkcc(f_LP, A_LP, b_LP, -Inf*ones(size(A_LP,2),1), Inf*ones(size(A_LP,2),1), repmat('U',size(A_LP,1),1), repmat('C',size(A_LP,2),1), 1, struct('savefilename','SimpleLP'));
%                 if status == 110 % infeasible LP
%                     D_x0 = Inf;
%                 end
%         end
%         
%     case 'CRN' % D_crn standardized discrepancy
%             
%         % Optimization problem
%         % D(x_0) = min_{(m, w) in P} n*(muhat - m)'*Sigma^{-1}*(muhat - m)
%          
%         % Formulate as quadratic program
% %         % Regularize small eigenvalues
% %         [V, D] = eig(sample_var);
% %         eig_diag = diag(D);
% %         eig_diag(eig_diag < 1e-3) = 1e-3;
% %         D = diag(eig_diag);
% %         sample_var = V*D*inv(V);
% %         sample_var = (sample_var + sample_var')/2; % Make symmetric
%         sample_var = sample_var + 0.01*speye(k); % Crude regularization
%         H_QP = [2*n_vec(1)*inv(sample_var(~zero_var_solns,~zero_var_solns)), zeros(k_bar,q); zeros(q,k_bar), zeros(q,q)];
%         f_QP = [-2*n_vec(1)*(sample_mean(~zero_var_solns)'/sample_var(~zero_var_solns,~zero_var_solns))'; zeros(q,1)];
%         A_QP = [A, C];
%         b_QP = b;
%         opt_val_offset = n_vec(1)*sample_mean(~zero_var_solns)'*(sample_var(~zero_var_solns,~zero_var_solns)\sample_mean(~zero_var_solns));
%         
% %         % Eigenvalue decomposition
% %         [V, D] = eig(sample_var);
% %         num_zeros = sum(diag(D) < 1e-6);
% %         Htilde = V(:,1 + num_zeros:end)*sqrt(D(1 + num_zeros:end, 1 + num_zeros:end));
% %         H_QP = [eye(k-num_zeros), zeros(k-num_zeros,q); zeros(q,k-num_zeros), zeros(q,q)];
% %         H_QP = H_QP + 0.00001*speye(k-num_zeros + q);
% %         f_QP = [-sqrt(2*n_vec(1))*((sample_mean)'*pinv(Htilde'))'; zeros(q, 1)];
% %         A_QP = 1/(sqrt(2*n_vec(1)))*[A, C]*[Htilde, zeros(size(Htilde,1), q); zeros(q, size(Htilde,2)), eye(q)];
% %         b_QP = b;
% %         opt_val_offset = n_vec(1)*sample_mean'*pinv(Htilde')*pinv(Htilde)*sample_mean;
% %         
%         
% %         % Schur decomposition approach
% %         [U, T] = schur(sample_var);
% %         H_QP = [2*n_vec(1)*inv(T), zeros(k,q); zeros(q,k), zeros(q,q)];
% %         f_QP = [-2*n_vec(1)*(sample_mean'*inv(U')*inv(T))'; zeros(q,1)];
% %         A_QP = [A, C]*U;
% %         b_QP = b;
% %         opt_val_offset = n_vec(1)*sample_mean'*inv(U')*inv(T)*inv(U)*sample_mean;
% %         
% %         % Moore-Penrose pseudoinverse        
% %         % Formulate as quadratic program
% %         H_QP = [2*n_vec(1)*pinv(sample_var(~zero_var_solns,~zero_var_solns)), zeros(k_bar,q); zeros(q,k_bar), zeros(q,q)];
% %         f_QP = [-2*n_vec(1)*(pinv(sample_var(~zero_var_solns,~zero_var_solns))*sample_mean(~zero_var_solns)); zeros(q,1)];
% %         A_QP = [A, C];
% %         b_QP = b;
% %         opt_val_offset = n_vec(1)*sample_mean(~zero_var_solns)'*pinv(sample_var(~zero_var_solns,~zero_var_solns))*sample_mean(~zero_var_solns);
% 
%         % Solve quadratic program (suppress outputs) %'Display','none',
%         options = optimoptions('quadprog','Display','none','MaxIter',500); %,'OptimalityTolerance',10^(-3)); % Default tolerance = 1e-8
%         [~, f_val, exitflag] = quadprog(H_QP, f_QP, A_QP, b_QP, [], [], [], [], [], options);
%         D_x0 = f_val + opt_val_offset;
% %         if exitflag ~= 1
% %             fprintf('\nExit flag = %d. Outputting %f.',exitflag, D_x0)
% %         end
%             
%         if exitflag == -2 % infeasible QP
%             D_x0 = Inf;
%         end
%             
% end % end switch

end