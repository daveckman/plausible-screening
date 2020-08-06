%load('tandem_M=1_iid_ellinf_convex.mat')
%load('tandem_M=1_iid_ell1_convex_budget50.mat')
load('tandem_M=1_iid_ellinf_convex_budget50.mat')
%load('tandem_M=1_iid_ell2_convex.mat')

[sort_D_x0s, rerank] = sort(D_x0s);

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)
xlim([0, card_feas_region])
%ylim([0,200])
xlabel('Sorted Solution Index','interpreter','latex')
ylabel('$D(x_0, \widehat{\mu}, \widehat{\Sigma}, n)$','interpreter','latex')
%title(plt_title)
hold on
plot(1:card_feas_region, sort_D_x0s, 'b-', 'LineWidth', 2);
line([0, card_feas_region], [D_cutoffs(discrep_index), D_cutoffs(discrep_index)], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5)
hold off

%print(['sorted_min_discrep_tandem_ellinf'],'-dpng','-r500')

%%
% Three Colors
load('tandem_M=1_iid_ellinf_convex.mat')

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
%[~, z, exit_flag] = linprog(f_LP, A_LP, b_LP, R(1:end-1,:), beq_LP(1:end-1), [], [], []);

%% Check if any 5 points are co-planar
%combos = combnk(1:100,5); % too big to enumerate
combos = combnk(1:50,5);
for i = 1:size(combos,1)
    vectors_to_test = exp_set(combos(i,:),:);
    vector_diffs = vectors_to_test - vectors_to_test(1,:);
    rank_check = rank(vector_diffs(2:5,:));
    if rank_check <= 2
        disp(vectors_to_test)
        break
    end
end

histogram(exp_set(:,1), 25);

%%

rand_exp_set = feas_region(randi([0, size(feas_region, 1)], 100, 1),:);
