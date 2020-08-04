function SX_indicators = construct_det_subset(feas_region, exp_set, true_mean, fn_props, prop_params)
% Construct the deterministic subset S(X) for known objective function
% values.

PO_info_handle = str2func(strcat('setup_',fn_props));
card_feas_region = size(feas_region, 1);
SX_indicators = zeros(card_feas_region, 1);

[~, exp_set_indices] = ismember(exp_set, feas_region, 'rows'); 
true_mean_exp_set = true_mean(exp_set_indices,:);

for l = 1:card_feas_region
    
    x0 = feas_region(l,:);
    
    [A, C, b] = PO_info_handle(x0, exp_set, prop_params);

    p = size(A,1); % Number of constraints
    k = size(A,2); % Number of solutions in experimental set
    q = size(C,2); % Number of unprojected components
    
    if isempty(C) % if explicit form of M(x0) is known, skip the optimization
    
        if sum(A*true_mean_exp_set <= b) == p % if all constraints are satisfied
            z = 0; % include x0 in S^poly
        else
            z = -Inf;
        end
    
    else % Formulate and solve a linear program (via Farkas' lemma)

        f_LP = b - A*true_mean_exp_set;
        A_LP = -speye(p);
        b_LP = zeros(p,1);
        Aeq_LP = C';
        beq_LP = zeros(q,1);

        %switch LP_solver_string
        %    case 'MATLAB'
        options = optimoptions('linprog','Display','none');
        [~, z, exit_flag] = linprog(f_LP, A_LP, b_LP, Aeq_LP, beq_LP, [], [], [], options);

        %Check if linear program was unbounded.
        if exit_flag == -3
            z = -Inf;
        end

%             case 'glpk'
%                 [~, z, exit_flag] = glpk(f_LP, [A_LP; Aeq_LP], [b_LP; beq_LP], -Inf*ones(size(A_LP,2),1), Inf*ones(size(A_LP,2),1), [repmat('U',size(A_LP,1),1); repmat('S',size(Aeq_LP,1),1)], repmat('C',size(A_LP,2),1), 1, struct('savefilename','SimpleLP'));
% 
%                 %Check if linear program was unbounded.
%                 if exit_flag == 6 || exit_flag == 111
%                     z = -Inf;
%                 end
%         end

    end

    SX_indicators(l) = (z >= 0); 
end

