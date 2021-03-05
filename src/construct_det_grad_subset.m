function SGX_indicators = construct_det_grad_subset(feas_region, exp_set, true_mean, true_grad) %, fn_props, prop_params)
% Construct the deterministic subset SG(X) for known objective function
% and gradient values.

% true_grad is a (|feasregion| x d) matrix, i.e., gradients are row vectors

% FOR CONVEX CASE

%PO_info_handle = str2func(strcat('setup_',fn_props));
k = size(exp_set, 1);
card_feas_region = size(feas_region, 1);
SGX_indicators = zeros(card_feas_region, 1);

[~, exp_set_indices] = ismember(exp_set, feas_region, 'rows'); 
true_mean_exp_set = true_mean(exp_set_indices,:);
true_grad_exp_set = true_grad(exp_set_indices,:); 

%
% SG(X) = {x0 in feasregion:
%          s.t. max_i {mu(x_i) - (x_i - x_0)^T*grad(x_i)} <= min_i mu(x_i)} 

for l = 1:card_feas_region
    
    x0 = feas_region(l,:);
    
    % WORKS FOR 1-DIM. STILL NEED TO FIX
    leftmax = max(true_mean_exp_set - (exp_set - repmat(x0, k, 1)).*true_grad_exp_set);
    rightmin = min(true_mean_exp_set);
    
    SGX_indicators(l) = (leftmax < rightmin);
    
%     
%     [A, C, b] = PO_info_handle(x0, exp_set, prop_params);
% 
%     p = size(A,1); % Number of constraints
%     k = size(A,2); % Number of solutions in experimental set
%     q = size(C,2); % Number of unprojected components
%     
%     if isempty(C) % if explicit form of M(x0) is known, skip the optimization
%     
%         if sum(A*true_mean_exp_set <= b) == p % if all constraints are satisfied
%             z = 0; % include x0 in S^poly
%         else
%             z = -Inf;
%         end
%     
%     else % Formulate and solve a linear program (via Farkas' lemma)
% 
%         f_LP = b - A*true_mean_exp_set;
%         A_LP = -speye(p);
%         b_LP = zeros(p,1);
%         Aeq_LP = C';
%         beq_LP = zeros(q,1);
% 
%         %switch LP_solver_string
%         %    case 'MATLAB'
%         options = optimoptions('linprog','Display','none');
%         [~, z, exit_flag] = linprog(f_LP, A_LP, b_LP, Aeq_LP, beq_LP, [], [], [], options);
% 
%         %Check if linear program was unbounded.
%         if exit_flag == -3
%             z = -Inf;
%         end
% 
%     end

%    SX_indicators(l) = (z >= 0); 

end

