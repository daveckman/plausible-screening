% When calculating S^poly, linprog sometimes "loses" feasibility.

clear
clc

addpath('../../tools')
addpath('../../src')
%addpath('../tools')
%addpath('../src')

%crunch_cluster = parcluster;
%maxNumCompThreads(4);

%load('tandem_M=1_iid_ellinf_mode=2_convex.mat')
%D_cutoff = D_cutoffs(3);

load('tandem_M=1_iid_ell2_mode=2_convex.mat')
D_cutoff = D_cutoffs(2);

bug_soln_indices = find(zs == Inf);

temp_zs = zeros(length(bug_soln_indices),1);

solns_to_fix = feas_region(bug_soln_indices,:);

parfor i = 1:length(bug_soln_indices)
%parfor (i = 1:length(bug_soln_indices), crunch_cluster)
    
    if mod(i,100) == 0
        fprintf('Solution %d.\n', i)
    end
    
        x0 = solns_to_fix(i,:)
        
        addpath('..\..\src')
        [A, C, b] = setup_convex(x0, exp_set, prop_params);

        p = size(A,1); % Number of constraints
        k = size(A,2); % Number of solutions in experimental set
        q = size(C,2); % Number of unprojected components

        % d infinity
        %bprime = b + D_cutoff * abs(A)*sqrt(sample_var./n_vec);

        % d 2
        bprime = b + sqrt(D_cutoff * (A.^2)*(sample_var./n_vec)); 
        
        f_LP = bprime - A*sample_mean;
        A_LP = -speye(p);
        b_LP = zeros(p,1);
        Aeq_LP = C';
        beq_LP = zeros(q,1);

        %[~, z, exit_flag] = linprog(f_LP, A_LP, b_LP, Aeq_LP, beq_LP, [], [], []);

        [~, z, exit_flag] = glpk(f_LP, [A_LP; Aeq_LP], [b_LP; beq_LP], -Inf*ones(size(A_LP,2),1), Inf*ones(size(A_LP,2),1), [repmat('U',size(A_LP,1),1); repmat('S',size(Aeq_LP,1),1)], repmat('C',size(A_LP,2),1), 1, struct('savefilename','SimpleLP'));

        %Check if linear program was unbounded.
        if exit_flag == 6 || exit_flag == 111
            z = -Inf;
        end

        temp_zs(i) = z;
end

zs(bug_soln_indices) = temp_zs;
S_poly_indicators(bug_soln_indices) = (temp_zs >= 0);

%rmpath('../tools')
%rmpath('../src')
rmpath('../../tools')
rmpath('../../src')

%save('tandem_M=1_iid_ellinf_mode=2_convex_amended.mat')
save('tandem_M=1_iid_ell2_mode=2_convex_amended.mat')