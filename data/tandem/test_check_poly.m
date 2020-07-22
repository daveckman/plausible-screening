load('tandem_M=1_iid_convex_budget_25.mat')

[sort_D_x0s, rerank] = sort(D_x0s);

group = zeros(size(D_x0s));
group(zs(rerank) == -Inf) = 1;
group(zs(rerank) == 0) = 2;
group(zs(rerank) == Inf) = 3;
groupColors = [0 0 1; 0 0 0; 1 0 0]; 
gscatter(1:size(D_x0s,1), sort_D_x0s, group, groupColors,'.')

xlabel('Solutions')
ylabel('Minimum Standardized Discrepancy')
title('Sorted Minimum Standardized Discrepancy')
line([0, card_feas_region], [D_cutoff_dinf, D_cutoff_dinf])
legend off
%intersect(feas_region(S_indicators_dinf==1,:), feas_region(S_poly_indicators_dinf==1,:),'rows')

%% Reproduce PO bug

bug_soln_index = find(zs == Inf, 1);
x0 = feas_region(bug_soln_index,:);

addpath('..\..\src')
[A, C, b] = setup_convex(x0, exp_set, prop_params);

p = size(A,1); % Number of constraints
k = size(A,2); % Number of solutions in experimental set
q = size(C,2); % Number of unprojected components

bprime = b + D_cutoff_dinf * abs(A)*sqrt(sample_var./n_vec);

f_LP = bprime - A*sample_mean;
A_LP = -speye(p);
b_LP = zeros(p,1);
Aeq_LP = C';
beq_LP = zeros(q,1);

[~, z, exit_flag] = linprog(f_LP, A_LP, b_LP, Aeq_LP, beq_LP, [], [], []);
R = rref(Aeq_LP);
[~, z, exit_flag] = linprog(f_LP, A_LP, b_LP, R(1:end-1,:), beq_LP(1:end-1), [], [], []);

