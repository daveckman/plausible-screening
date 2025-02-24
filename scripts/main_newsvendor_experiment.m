% Script for testing the interface

clear;
clc;

add_rm_paths('add');

problem_string = 'newsvendor';
[oracle_string, oracle_n_rngs, feas_region, exp_set, k, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string] = init_problem(problem_string);

check_exceptions(discrep_string, fn_props, n_vec)

card_feas_region = size(feas_region, 1);
per_soln_samplesize = sum(n_vec)/card_feas_region;

if floor(per_soln_samplesize) == per_soln_samplesize && per_soln_samplesize >= 2
    n_vec_SS = per_soln_samplesize*ones(card_feas_region, 1);
else
    fprintf('Total size of %d cannot be allocated equally across %d feasible solutions.\n', sum(n_vec), card_feas_region);
end

%% RUN MACROREPLICATIONS

M = 200; % Number of macroreplications

% Initialize data storage
S_indicators_d1 = zeros(card_feas_region, M);
S_indicators_d2 = zeros(card_feas_region, M);
S_indicators_dinf = zeros(card_feas_region, M);
%S_indicators_dcrn = zeros(card_feas_region, M);
S_poly_indicators_d1 = zeros(card_feas_region, M);
S_poly_indicators_d2 = zeros(card_feas_region, M);
S_poly_indicators_dinf = zeros(card_feas_region, M);
%S_poly_indicators_dcrn = zeros(card_feas_region, M);
SS_indicators = zeros(card_feas_region, M);

parfor m = 1:M
    
    fprintf('Running macrorep %d of %d.\n', m, M)
    
    % SAMPLING
    
%     % Generate data and calculate summary statistics
%     fprintf('Generating sample data for plausible optima...\n')
%     [sample_mean, sample_var] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, 'ell1');
% 
%     % SCREENING (USING DIFFERENT DISCREPANCIES)
%     discrep_string = 'ell1';
%     fprintf('Screening solutions for %s discrepancy...\n', discrep_string)
%     [S_indicators_d1(:,m), D_x0s, S_poly_indicators_d1(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string);
%     fprintf('\nResults of PO screening\n-------------------------------------------------------\n')
%     fprintf('\tstandardized discrepancy: \t\t\t\t%s\n', discrep_string)
%     fprintf('\t# of solutions in PO subset: \t\t\t%d\n', sum(S_indicators_d1(:,m)))
%     fprintf('\t# of solutions in PO relaxed subset: \t%d\n\n', sum(S_poly_indicators_d1(:,m)))
% 
%     discrep_string = 'ell2';
%     fprintf('Screening solutions for %s discrepancy...\n', discrep_string)
%     [S_indicators_d2(:,m), D_x0s, S_poly_indicators_d2(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string);
%     fprintf('\nResults of PO screening\n-------------------------------------------------------\n')
%     fprintf('\tstandardized discrepancy: \t\t\t\t%s\n', discrep_string)
%     fprintf('\t# of solutions in PO subset: \t\t\t%d\n', sum(S_indicators_d2(:,m)))  
%     fprintf('\t# of solutions in PO relaxed subset: \t%d\n\n', sum(S_poly_indicators_d2(:,m)))
% 
%     discrep_string = 'ellinf';
%     fprintf('Screening solutions for %s discrepancy...\n', discrep_string)
%     [S_indicators_dinf(:,m), D_x0s, S_poly_indicators_dinf(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string);
%     fprintf('\nResults of PO screening\n-------------------------------------------------------\n')
%     fprintf('\tstandardized discrepancy: \t\t\t\t%s\n', discrep_string)
%     fprintf('\t# of solutions in PO subset: \t\t\t%d\n', sum(S_indicators_dinf(:,m)))
%     fprintf('\t# of solutions in PO relaxed subset: \t%d\n\n', sum(S_poly_indicators_dinf(:,m)))

    % SCREENING (USING EXTENDED SCREEN TO THE BEST)

    % Generate data and calculate summary statistics
    fprintf('Generating sample data for subset selection...\n')
    [sample_mean_SS, sample_var_SS] = generate_data(m, oracle_string, oracle_n_rngs, feas_region, n_vec_SS, 'ell1');
    
    fprintf('Screening solutions for extended screen-to-the-best...\n')
    [SS_indicators(:,m)] = ExtSTB(card_feas_region, sample_mean_SS, sample_var_SS, n_vec_SS, alpha);
    fprintf('\nResults of Extended Screen-to-the-Best\n-------------------------------------------------------\n')
    fprintf('\t# of solutions in ExtSTB subset: \t\t\t%d\n\n', sum(SS_indicators(:,m)))
    
%     % Generate data and calculate summary statistics
%     fprintf('Generating sample data in parallel...\n')
%     [sample_mean, sample_var] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, 'CRN');
% 
%     discrep_string = 'CRN';
%     fprintf('Screening solutions in parallel for %s discrepancy...\n', discrep_string)
%     [S_indicators_dcrn(:,m), ~, S_poly_indicators_dcrn(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string);
%     fprintf('\nResults of PO screening\n-------------------------------------------------------\n')
%     fprintf('\tstandardized discrepancy: \t\t\t\t%s\n', discrep_string)
%     fprintf('\t# of solutions in PO subset: \t\t\t%d\n', sum(S_indicators_dcrn(:,m)))
%     fprintf('\t# of solutions in PO relaxed subset: \t%d\n\n', sum(S_poly_indicators_dcrn(:,m)))
%     
end

%%
% ESTIMATE TRUE OBJECTIVE FUNCTION
M_MC = 1000;
outputs = zeros(card_feas_region, M_MC);

oracle_rngs = {RandStream.create('mrg32k3a')};

for i = 1:card_feas_region
    % Extract solution x_i and sample size n_i
    x_i = feas_region(i,:);

    % Reset each stream to substream 1
    for r = 1:oracle_n_rngs
        oracle_rng = oracle_rngs{r};
        oracle_rng.Substream = 1;
    end

    % Take n_i replications at x_i
    outputs(i,:) = newsvendor_oracle(oracle_rngs, x_i, M_MC);
end

% Calculate summary statistics
est_true_mean = mean(outputs,2);

% Calculate optimal solution
cost = 3; % per unit order cost 
sell_price = 9; % per unit sale revenue
salvage = 1; % per unit salvage revenue
lambda = 50; % average daily demand
opt_sol = ceil(poissinv((sell_price - cost)/(sell_price - salvage), lambda));

%%
% PLOTTING SUBSETS

figure
subplot(1, 2, 1);

axis square
xlim([min(feas_region), max(feas_region)])
ylim([0,1])
xlabel('Order Quantity','interpreter','latex')
ylabel('$P(x_0 \in S)$','interpreter','latex')
title('$d^1$ discrepancy (standard PO)', 'interpreter', 'latex')

hold on
yyaxis left
C = sum(S_indicators_d1,2)'/M;
scatter(feas_region(:,1), C, 'bo','markerfacecolor','b')
plot(exp_set(:,1),zeros(1,k),'kd','markerfacecolor','k')
scatter(opt_sol, 0, 'gs', 'markerfacecolor', 'g')
yyaxis right
plot(feas_region(:,1),est_true_mean,'r-')
hold off

subplot(1, 2, 2);

axis square
xlim([min(feas_region), max(feas_region)])
ylim([0,1])
xlabel('Order Quantity','interpreter','latex')
ylabel('$P(x_0 \in S)$','interpreter','latex')
title('$d^1$ discrepancy (relaxed PO)', 'interpreter', 'latex')

hold on
yyaxis left
C = sum(S_poly_indicators_d1,2)'/M;
scatter(feas_region(:,1), C, 'bo','markerfacecolor','b')
plot(exp_set(:,1),zeros(1,k),'kd','markerfacecolor','k')
scatter(opt_sol, 0, 'gs', 'markerfacecolor', 'g')
yyaxis right
plot(feas_region(:,1),est_true_mean,'r-')
hold off

figure;
subplot(1, 2, 1);


axis square
xlim([min(feas_region), max(feas_region)])
ylim([0,1])
xlabel('Order Quantity','interpreter','latex')
ylabel('$P(x_0 \in S)$','interpreter','latex')
title('$d^2$ discrepancy (standard PO)', 'interpreter', 'latex')

hold on
yyaxis left
C = sum(S_indicators_d2,2)'/M;
scatter(feas_region(:,1), C, 'bo','markerfacecolor','b')
plot(exp_set(:,1),zeros(1,k),'kd','markerfacecolor','k')
scatter(opt_sol, 0, 'gs', 'markerfacecolor', 'g')
yyaxis right
plot(feas_region(:,1),est_true_mean,'r-')
hold off

subplot(1, 2, 2);

axis square
xlim([min(feas_region), max(feas_region)])
ylim([0,1])
xlabel('Order Quantity','interpreter','latex')
ylabel('$P(x_0 \in S)$','interpreter','latex')
title('$d^2$ discrepancy (relaxed PO)', 'interpreter', 'latex')

hold on
yyaxis left
C = sum(S_poly_indicators_d2,2)'/M;
scatter(feas_region(:,1), C, 'bo','markerfacecolor','b')
plot(exp_set(:,1),zeros(1,k),'kd','markerfacecolor','k')
scatter(opt_sol, 0, 'gs', 'markerfacecolor', 'g')
yyaxis right
plot(feas_region(:,1),est_true_mean,'r-')
hold off

figure;
subplot(1, 2, 1)

axis square
xlim([min(feas_region), max(feas_region)])
ylim([0,1])
xlabel('Order Quantity','interpreter','latex')
ylabel('$P(x_0 \in S)$','interpreter','latex')
title('$d^\infty$ discrepancy (standard PO)', 'interpreter', 'latex')

hold on
yyaxis left
C = sum(S_indicators_dinf,2)'/M;
scatter(feas_region(:,1), C, 'bo','markerfacecolor','b')
plot(exp_set(:,1),zeros(1,k),'kd','markerfacecolor','k')
scatter(opt_sol, 0, 'gs', 'markerfacecolor', 'g')
yyaxis right
plot(feas_region(:,1),est_true_mean,'r-')
hold off

subplot(1, 2, 2);

axis square
xlim([min(feas_region), max(feas_region)])
ylim([0,1])
xlabel('Order Quantity','interpreter','latex')
ylabel('$P(x_0 \in S)$','interpreter','latex')
title('$d^\infty$ discrepancy (relaxed PO)', 'interpreter', 'latex')

hold on
yyaxis left
C = sum(S_poly_indicators_dinf,2)'/M;
scatter(feas_region(:,1), C, 'bo','markerfacecolor','b');
plot(exp_set(:,1),zeros(1,k),'kd','markerfacecolor','k')
scatter(opt_sol, 0, 'gs', 'markerfacecolor', 'g')
yyaxis right
plot(feas_region(:,1),est_true_mean,'r-')
hold off

%% Subset selection

figure;

axis square
xlim([min(feas_region), max(feas_region)])
ylim([0,1])
xlabel('Order Quantity','interpreter','latex')
ylabel('$P(x_0 \in S)$','interpreter','latex')
title('Extended Screen-to-the-Best', 'interpreter', 'latex')

hold on
yyaxis left
C = sum(SS_indicators,2)'/M;
scatter(feas_region(:,1), C, 'bo','markerfacecolor','b')
plot(exp_set(:,1),zeros(1,k),'kd','markerfacecolor','k')
scatter(opt_sol, 0, 'gs', 'markerfacecolor', 'g')
yyaxis right
plot(feas_region(:,1),est_true_mean,'r-')
hold off

% figure;
% subplot(1, 2, 1)
% 
% axis square
% xlim([min(feas_region), max(feas_region)])
% ylim([0,1])
% xlabel('Order Quantity','interpreter','latex')
% ylabel('$P(x_0 \in S)$','interpreter','latex')
% title('$d^{CRN}$ discrepancy (standard PO)', 'interpreter', 'latex')
% 
% hold on
% yyaxis left
% C = sum(S_indicators_dcrn,2)'/M;
% scatter(feas_region(:,1), C, 'bo','markerfacecolor','b')
% plot(exp_set(:,1),zeros(1,k),'kd','markerfacecolor','k')
% scatter(opt_sol, 0, 'gs', 'markerfacecolor', 'g')
% yyaxis right
% plot(feas_region(:,1),est_true_mean,'r-')
% hold off
% 
% subplot(1, 2, 2);
% 
% axis square
% xlim([min(feas_region), max(feas_region)])
% ylim([0,1])
% xlabel('Order Quantity','interpreter','latex')
% ylabel('$P(x_0 \in S)$','interpreter','latex')
% title('$d^{CRN}$ discrepancy (relaxed PO)', 'interpreter', 'latex')
% 
% hold on
% yyaxis left
% C = sum(S_poly_indicators_dcrn,2)'/M;
% scatter(feas_region(:,1), C, 'bo','markerfacecolor','b');
% plot(exp_set(:,1),zeros(1,k),'kd','markerfacecolor','k')
% scatter(opt_sol, 0, 'gs', 'markerfacecolor', 'g')
% yyaxis right
% plot(feas_region(:,1),est_true_mean,'r-')
% hold off

%% PLOTTING SUBSET SIZES (ECDFS)

% Compute subset sizes
subsize_d1 = sum(S_indicators_d1, 1);
subsize_d1_poly = sum(S_poly_indicators_d1, 1);
subsize_d2 = sum(S_indicators_d2, 1);
subsize_d2_poly = sum(S_poly_indicators_d2, 1);
subsize_dinf = sum(S_indicators_dinf, 1);
subsize_dinf_poly = sum(S_poly_indicators_dinf, 1);
subsize_SS = sum(SS_indicators, 1);

% Compute empirical cdf of expected sample size
ecdf_d1 = zeros(1, card_feas_region);
ecdf_d1_poly = zeros(1, card_feas_region);
ecdf_d2 = zeros(1, card_feas_region);
ecdf_d2_poly = zeros(1, card_feas_region);
ecdf_dinf = zeros(1, card_feas_region);
ecdf_dinf_poly = zeros(1, card_feas_region);
ecdf_SS = zeros(1, card_feas_region);

for j = 1:card_feas_region
    ecdf_d1(j) = mean(subsize_d1 <= j);
    ecdf_d1_poly(j) = mean(subsize_d1_poly <= j);
    ecdf_d2(j) = mean(subsize_d2 <= j);
    ecdf_d2_poly(j) = mean(subsize_d2_poly <= j);
    ecdf_dinf(j) = mean(subsize_dinf <= j);
    ecdf_dinf_poly(j) = mean(subsize_dinf_poly <= j);
    ecdf_SS(j) = mean(subsize_SS <= j);
end

figure
hold on

stairs([1:card_feas_region], ecdf_d1, 'b-', 'LineWidth', 2);
stairs([1:card_feas_region], ecdf_d1_poly, 'b:','LineWidth', 2);
stairs([1:card_feas_region], ecdf_d2, 'r-', 'LineWidth', 2);
stairs([1:card_feas_region], ecdf_d2_poly, 'r:', 'LineWidth', 2);
stairs([1:card_feas_region], ecdf_dinf, 'g-', 'LineWidth', 2);
stairs([1:card_feas_region], ecdf_dinf_poly, 'g:', 'LineWidth', 2);
stairs([1:card_feas_region], ecdf_SS, 'c-', 'LineWidth', 2);

xlim([0,card_feas_region])

xlabel('Subset Size ($s$)', 'Interpreter', 'latex')
ylabel('P($|S| \leq s$)', 'Interpreter', 'latex')
leg1 = legend('$d^1$ standard', '$d^1$ relaxed', '$d^2$ standard', '$d^2$ relaxed', '$d^{\infty}$ standard', '$d^{\infty}$ relaxed', 'ExtSTB', 'location', 'northwest');
legend boxoff
set(leg1,'Interpreter','latex');
set(gca, 'FontSize', 14, 'LineWidth', 2)
set(gcf,'units','inches','position',[1,1,9,4.5])

hold off

add_rm_paths('remove');
