function z = check_poly_feas(discrep_string, A, C, b, sample_mean, sample_var, n_vec, D_cutoff, LP_solver_string)

% 1. Take polyhedral representation of P.
% 2. Offset the rhs vector to account for the uncertainty about the true function.
% 3. Check whether sample_mean belongs to the relaxed polyhedron by either
%    directly checking the constraints of solving a linear program (from 
%    Farkas' lemma).

% Determine dimensions of polyhedral representation of M(x_0) = proj_m(P)
% where P = {(m, w): A*m + C*w <= b}
p = size(A,1); % Number of constraints
k = size(A,2); % Number of solutions in experimental set
q = size(C,2); % Number of unprojected components

% linprog(f_LP, A_LP, b_LP, Aeq_LP, beq_LP) solves a linear program of the form
%   min_x f_LP'*x s.t. A_LP*x <= b_LP and Aeq_LP*x = beq_LP where
%       f_LP is a column vector
%       A_LP is a matrix
%       b_LP is a column vector

% The linear program for checking feasibility is
%   min_y (bprime-A*sample_mean)^T*y s.t. y >= 0 and C^T*y = 0

bprime = zeros(p,1);

switch discrep_string
    case 'ell1' % D_1 standardized discrepancy
        
        bprime = b + D_cutoff * max(sqrt(sample_var./n_vec)'.*abs(A), [], 2);
        
    case 'ell2' % D_2 standardized discrepancy
       
        bprime = b + sqrt(D_cutoff * (A.^2)*(sample_var./n_vec)); 
        
    case 'ellinf' % D_inf standardized discrepancy

        bprime = b + D_cutoff * abs(A)*sqrt(sample_var./n_vec);
        
    case 'CRN' % D_crn standardized discrepancy
        
        term = 1/n_vec(1)*sum((A*sample_var).*A,2);
        term(term < 0) = 0; % avoid negative signed zeros leading to imag #s
        bprime = b + sqrt(D_cutoff * term);
        
end % end switch

if isempty(C) % if explicit form of M(x0) is known, skip the optimization
    
    if sum(A*sample_mean <= bprime) == p % if all constraints are satisfied
        z = 0; % include x0 in S^poly
    else
        z = -Inf;
    end
    
else % Formulate and solve a linear program 
    
    % extra variable approach    
    % Solve z = max_{w, rho} rho s.t. Cw + rho*1_p <= b' - A*muhat
    % flip objective because of maximization
    f_LP = [zeros(q,1);-1];
    A_LP = [C, ones(p, 1)];
    b_LP = bprime - A*sample_mean;

    options = optimoptions('linprog','Display','none');
    [~, z, exit_flag] = linprog(f_LP, A_LP, b_LP,[], [], [], [], [], options);
        
    %Check if linear program was unbounded.
    if exit_flag == -3
        z = -Inf;
    end
    
    % Flip sign back
    z = -z;
    
% Farkas' Lemma approach
%     f_LP = bprime - A*sample_mean;
%     A_LP = -speye(p);
%     b_LP = zeros(p,1);
%     Aeq_LP = C';
%     beq_LP = zeros(q,1);
% 
%     switch LP_solver_string
%         case 'MATLAB'
%             options = optimoptions('linprog','Display','none');
%             [~, z, exit_flag] = linprog(f_LP, A_LP, b_LP, Aeq_LP, beq_LP, [], [], [], options);
%             
%             %Check if linear program was unbounded.
%             if exit_flag == -3
%                 z = -Inf;
%             end
%             
%         case 'glpk'
%             [~, z, exit_flag] = glpk(f_LP, [A_LP; Aeq_LP], [b_LP; beq_LP], -Inf*ones(size(A_LP,2),1), Inf*ones(size(A_LP,2),1), [repmat('U',size(A_LP,1),1); repmat('S',size(Aeq_LP,1),1)], repmat('C',size(A_LP,2),1), 1, struct('savefilename','SimpleLP'));
%             
%             %Check if linear program was unbounded.
%             if exit_flag == 6 || exit_flag == 111
%                 z = -Inf;
%             end
%     end
    
end
    
end