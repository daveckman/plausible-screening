% Script for testing the interface

clear;
clc;

add_rm_paths('add');

problem_string = 'synthetic';
[oracle_string, oracle_n_rngs, feas_region, exp_set, k, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string] = init_problem(problem_string);

check_exceptions(discrep_string, fn_props, n_vec)

card_feas_region = size(feas_region, 1);
%per_soln_samplesize = sum(n_vec)/card_feas_region;

% if floor(per_soln_samplesize) == per_soln_samplesize && per_soln_samplesize >= 2
%     n_vec_SS = per_soln_samplesize*ones(card_feas_region, 1);
% else
%     fprintf('Total size of %d cannot be allocated equally across %d feasible solutions.\n', sum(n_vec), card_feas_region);
% end

%% RUN MACROREPLICATIONS

M = 100; % Number of macroreplications

% Calculate cutoffs for PS
D_cutoff_d2 = calc_cutoff(k, n_vec, alpha, 'ell2');
d = size(feas_region, 2);
D_cutoff_grad = calc_grad_cutoff(k, d, n_vec, alpha);
D_cutoff_gradinf = calc_gradinf_cutoff(k, d, n_vec, alpha);

% Initialize data storage
S_indicators_d2 = zeros(card_feas_region, M);
D_x0s = zeros(card_feas_region, M);
S_poly_indicators_d2 = zeros(card_feas_region, M);

S_indicators_grad = zeros(card_feas_region, M);
D_x0s_grad = zeros(card_feas_region, M);
S_poly_indicators_grad = zeros(card_feas_region, M);

Sonlygrad_d2_indicators = zeros(card_feas_region, M);
Sonlygrad_dinf_indicators = zeros(card_feas_region, M);

print_problem_header(problem_string, feas_region, exp_set, fn_props)

parfor m = 1:M
    
    fprintf('Running macrorep %d of %d.\n', m, M)
    
    % SAMPLING
    
    % Generate data and calculate summary statistics
    fprintf('Generating sample data for plausible optima...\n')
    [sample_mean, sample_var, sample_mean_grad, sample_full_cov] = generate_grad_data(m, 'synthetic_grad_oracle', oracle_n_rngs, exp_set, n_vec, 'grad');

    %[sample_mean, sample_var] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, 'ell1');

%     % Original Plausible Screening
%     discrep_string = 'ell2';
%     fprintf('Screening solutions for %s discrepancy...\n', discrep_string)
%     [S_indicators_d2(:,m), D_x0s(:,m), S_poly_indicators_d2(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ell2', D_cutoff_d2, fn_props, prop_params, LP_solver_string);
%     fprintf('\nResults of PS screening\n-------------------------------------------------------\n')
%     fprintf('\tstandardized discrepancy: \t\t\t\t%s\n', discrep_string)
%     fprintf('\t# of solutions in PS subset: \t\t\t%d\n', sum(S_indicators_d2(:,m)))  
%     fprintf('\t# of solutions in PS relaxed subset: \t%d\n\n', sum(S_poly_indicators_d2(:,m)))
% 
%     % Gradient Plausible Screening w/ Dgrad
%     [S_indicators_grad(:,m), D_x0s_grad(:,m), S_poly_indicators_grad(:,m)] = GPS_screen(feas_region, exp_set, sample_mean, sample_mean_grad, sample_full_cov, n_vec, D_cutoff_grad);
%     fprintf('\nResults of GPS screening\n-------------------------------------------------------\n')
%     fprintf('\tstandardized discrepancy: \t\t\t\t Gradients\n')
%     fprintf('\t# of solutions in GPS subset: \t\t\t%d\n', sum(S_indicators_grad(:,m)))  
%     fprintf('\t# of solutions in GPS relaxed subset: \t%d\n\n', sum(S_poly_indicators_grad(:,m)))
% 
%     % Gradient Plausible Screening w/ Dgradinf
%     [S_indicators_gradinf(:,m), D_x0s_gradinf(:,m), S_poly_indicators_gradinf(:,m)] = GPSinf_screen(feas_region, exp_set, sample_mean, sample_mean_grad, sample_full_cov, n_vec, D_cutoff_gradinf);
%     fprintf('\nResults of GPS screening (dinf) \n-------------------------------------------------------\n')
%     fprintf('\tstandardized discrepancy: \t\t\t\t Gradients\n')
%     fprintf('\t# of solutions in GPS inf subset: \t\t\t%d\n', sum(S_indicators_gradinf(:,m)))  
%     fprintf('\t# of solutions in GPS inf relaxed subset: \t%d\n\n', sum(S_poly_indicators_gradinf(:,m)))
%     
    % Relaxed Gradient Plausible Screening w/ only gradients
    [Sonlygrad_d2_indicators(:,m), Sonlygrad_dinf_indicators(:,m)] = RGPS_onlygrad_screen(feas_region, exp_set, sample_mean, sample_mean_grad, sample_full_cov, n_vec, D_cutoff_grad, D_cutoff_gradinf);
    fprintf('\nResults of RGPS screening (only gradients) \n-------------------------------------------------------\n')
    fprintf('\t# of solutions in RGPS (only gradients) d2 subset: \t\t\t%d\n', sum(Sonlygrad_d2_indicators(:,m)))  
    fprintf('\t# of solutions in RGPS (only gradients) dinf subset: \t%d\n\n', sum(Sonlygrad_dinf_indicators(:,m)))
end

%%
% CALCULATE TRUE OBJECTIVE FUNCTIONS AND TRUE GRADIENTS

% mu(x) = x^2 - 0.5*x + 1
% GRAD mu(x) = 2*x - 0.5
true_mean = feas_region.^2 - 0.5*feas_region + 1;
true_grad = 2*feas_region - 0.5;

%%
% COMPUTE S(X)
SX_indicators = construct_det_subset(feas_region, exp_set, true_mean, 'convex', '');

switch_on = find(diff(SX_indicators) == 1)+1;
switch_off = find(diff(SX_indicators) == -1);

if min(switch_off) < min(switch_on)
    switch_on = [1; switch_on];
end
if max(switch_off) < max(switch_on)
    switch_off = [switch_off; 200];
end

%% Plot single 1-D P(x0 in S) for d2

alpha = 0.05;
%load(['ctsnews_N=',num2str(N),'_K=',num2str(K),'_M=',num2str(M),'_iid_',fn_property,'.mat'],'exp_set');

grey_rgb = (192/256)*ones(1,3);
dark_grey_rgb = (128/256)*ones(1,3);
figure
set(gca, 'FontSize', 14, 'LineWidth', 2)

xlim([0,1])
ylim([0,1.005])
xlabel('Solution ($x_0$)','interpreter','latex')
ylabel('Probability of $x_0$ in Subset','interpreter','latex')
%title(string_names{1},'interpreter','latex')

hold on
plot(feas_region, (true_mean - min(true_mean))./(max(true_mean) - min(true_mean)), 'color', dark_grey_rgb, 'LineWidth', 1)
line([0,1], [1-alpha, 1-alpha], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5)
C3 = mean(S_indicators_d2,2); %sum(all_S_indicators_d2(:,:,1),2)'/M;
C5 = mean(S_poly_indicators_d2,2); 
%scatter(feas_region(:,1), C, 'ko','markerfacecolor','k')
h1 = plot(feas_region(:,1), C3, 'b-', 'LineWidth', 2);
h2 = plot(feas_region(:,1), C5, 'r-', 'LineWidth', 2);
plot(exp_set,zeros(1,size(exp_set,2)),'kx','markerfacecolor','k', 'MarkerSize', 12, 'LineWidth', 1)

for i = 1:length(switch_on)
    x = [switch_on(i)-.5, switch_off(i)+.5, switch_off(i)+.5, switch_on(i)-.5]/card_feas_region;
    y = [0, 0, 1, 1]; 
    p=patch(x,y,'b','LineStyle','none','FaceAlpha',0.2);
end

legend([h1, h2], 'PS w/ d2', 'RPS w/ d2', 'Location', 'southeast')

box on
hold off

%print(['inc_probs_ctsnews_',fn_property,'_d2_K=5_N=400'],'-dpng','-r500')

%%
% COMPUTE SG(X)

SGX_indicators = construct_det_grad_subset(feas_region, exp_set, true_mean, true_grad);

switch_on = find(diff(SGX_indicators) == 1)+1;
switch_off = find(diff(SGX_indicators) == -1);

if min(switch_off) < min(switch_on)
    switch_on = [1; switch_on];
end
if max(switch_off) < max(switch_on)
    switch_off = [switch_off; 200];
end

%% Plot single 1-D P(x0 in S) for dgrad

%K = 5;
%N = 400;
%M = 200;

alpha = 0.05;
%load(['ctsnews_N=',num2str(N),'_K=',num2str(K),'_M=',num2str(M),'_iid_',fn_property,'.mat'],'exp_set');

grey_rgb = (192/256)*ones(1,3);
dark_grey_rgb = (128/256)*ones(1,3);
figure
set(gca, 'FontSize', 14, 'LineWidth', 2)

xlim([0,1])
ylim([0,1.005])
xlabel('Solution ($x_0$)','interpreter','latex')
ylabel('Probability of $x_0$ in Subset','interpreter','latex')
%title(string_names{1},'interpreter','latex')

hold on
plot(feas_region, (true_mean - min(true_mean))./(max(true_mean) - min(true_mean)), 'color', dark_grey_rgb, 'LineWidth', 1)
line([0,1], [1-alpha, 1-alpha], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5)
C8 = mean(S_indicators_grad,2); %sum(all_S_indicators_d2(:,:,1),2)'/M;
C4 = mean(S_poly_indicators_grad,2); %sum(all_S_indicators_d2(:,:,1),2)'/M;
%scatter(feas_region(:,1), C, 'ko','markerfacecolor','k')
h4 = plot(feas_region(:,1), C4, 'r-', 'LineWidth', 2);
h3 = plot(feas_region(:,1), C8, 'b-', 'LineWidth', 2);
plot(exp_set,zeros(1,size(exp_set,2)),'kx','markerfacecolor','k', 'MarkerSize', 12, 'LineWidth', 1)

for i = 1:length(switch_on)
    x = [switch_on(i)-.5, switch_off(i)+.5, switch_off(i)+.5, switch_on(i)-.5]/card_feas_region;
    y = [0, 0, 1, 1]; 
    p=patch(x,y,'b','LineStyle','none','FaceAlpha',0.2);
end

legend([h3, h4], 'GPS w/ dgrad', 'RGPS w/ dgrad', 'Location', 'southeast')

box on
hold off

%% Plot single 1-D P(x0 in S) for dgradinf

%K = 5;
%N = 400;
%M = 200;

alpha = 0.05;
%load(['ctsnews_N=',num2str(N),'_K=',num2str(K),'_M=',num2str(M),'_iid_',fn_property,'.mat'],'exp_set');

grey_rgb = (192/256)*ones(1,3);
dark_grey_rgb = (128/256)*ones(1,3);
figure
set(gca, 'FontSize', 14, 'LineWidth', 2)

xlim([0,1])
ylim([0,1.005])
xlabel('Solution ($x_0$)','interpreter','latex')
ylabel('Probability of $x_0$ in Subset','interpreter','latex')
%title(string_names{1},'interpreter','latex')

hold on
plot(feas_region, (true_mean - min(true_mean))./(max(true_mean) - min(true_mean)), 'color', dark_grey_rgb, 'LineWidth', 1)
line([0,1], [1-alpha, 1-alpha], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5)
C9 = mean(S_indicators_gradinf,2); %sum(all_S_indicators_d2(:,:,1),2)'/M;
C10 = mean(S_poly_indicators_gradinf,2); %sum(all_S_indicators_d2(:,:,1),2)'/M;
%scatter(feas_region(:,1), C, 'ko','markerfacecolor','k')
h4 = plot(feas_region(:,1), C10, 'r-', 'LineWidth', 2);
h3 = plot(feas_region(:,1), C9, 'b-', 'LineWidth', 2);
plot(exp_set,zeros(1,size(exp_set,2)),'kx','markerfacecolor','k', 'MarkerSize', 12, 'LineWidth', 1)

for i = 1:length(switch_on)
    x = [switch_on(i)-.5, switch_off(i)+.5, switch_off(i)+.5, switch_on(i)-.5]/card_feas_region;
    y = [0, 0, 1, 1]; 
    p=patch(x,y,'b','LineStyle','none','FaceAlpha',0.2);
end

legend([h3, h4], 'GPS w/ dgradinf', 'RGPS w/ dgradinf', 'Location', 'southeast')

box on
hold off

%% Plot single 1-D P(x0 in S) for RGPS w/ only gradients (relaxation of M(x0))

alpha = 0.05;
%load(['ctsnews_N=',num2str(N),'_K=',num2str(K),'_M=',num2str(M),'_iid_',fn_property,'.mat'],'exp_set');

grey_rgb = (192/256)*ones(1,3);
dark_grey_rgb = (128/256)*ones(1,3);
figure
set(gca, 'FontSize', 14, 'LineWidth', 2)

xlim([0,1])
ylim([0,1.005])
xlabel('Solution ($x_0$)','interpreter','latex')
ylabel('Probability of $x_0$ in Subset','interpreter','latex')
%title(string_names{1},'interpreter','latex')

hold on
plot(feas_region, (true_mean - min(true_mean))./(max(true_mean) - min(true_mean)), 'color', dark_grey_rgb, 'LineWidth', 1)
line([0,1], [1-alpha, 1-alpha], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5)
C11 = mean(Sonlygrad_d2_indicators,2); %sum(all_S_indicators_d2(:,:,1),2)'/M;
C12 = mean(Sonlygrad_dinf_indicators,2); %sum(all_S_indicators_d2(:,:,1),2)'/M;
%scatter(feas_region(:,1), C, 'ko','markerfacecolor','k')
h4 = plot(feas_region(:,1), C12, 'r-', 'LineWidth', 2);
h3 = plot(feas_region(:,1), C11, 'b-', 'LineWidth', 2);
plot(exp_set,zeros(1,size(exp_set,2)),'kx','markerfacecolor','k', 'MarkerSize', 12, 'LineWidth', 1)

for i = 1:length(switch_on)
    x = [switch_on(i)-.5, switch_off(i)+.5, switch_off(i)+.5, switch_on(i)-.5]/card_feas_region;
    y = [0, 0, 1, 1]; 
    p=patch(x,y,'b','LineStyle','none','FaceAlpha',0.2);
end

legend([h3, h4], 'RGPS w/ d2 (only grad)', 'RGPS w/ dinf (only grad)', 'Location', 'southeast')

box on
hold off

%% THE REST IS NOT USED

%% 
% CHECK RGPS WITH PLUGGING IN TRUE GRADIENTS
%M = 100;

exact_grads = true_grad(exp_set,:);

D_grad = calc_grad_cutoff(k, d, n_vec, alpha);
D_d2 = calc_cutoff(k, n_vec, alpha, 'ell2');

S_test_indicators_grad = zeros(card_feas_region, M);

parfor m = 1:M
    
    fprintf('Running macrorep %d of %d.\n', m, M)
    
    % SAMPLING
    
    % Generate data and calculate summary statistics
    fprintf('Generating sample data for plausible optima...\n')
    [sample_mean, sample_var, sample_mean_grad, sample_full_cov] = generate_grad_data(m, 'cts_newsvendor_grad_oracle', oracle_n_rngs, exp_set, n_vec, 'grad');

    %[sample_mean, sample_var] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, 'ell1');
    S_testg_indicators_grad(:,m) = RGPSexact_screen(feas_region, exp_set, sample_mean, sample_var, exact_grads, n_vec, D_grad);
    S_test2_indicators_grad(:,m) = RGPSexact_screen(feas_region, exp_set, sample_mean, sample_var, exact_grads, n_vec, D_d2);
    fprintf('\nResults of RGPS screening (using exact gradients)\n-------------------------------------------------------\n')
    fprintf('\tstandardized discrepancy: \t\t\t\t Gradients\n')
    fprintf('\t# of solutions in GPS relaxed subset: \t%d\n\n', sum(S_testg_indicators_grad(:,m)))
    fprintf('\t# of solutions in GPS relaxed subset (d2 cutoff): \t%d\n\n', sum(S_test2_indicators_grad(:,m)))

end

%%

figure

alpha = 0.05;
%load(['ctsnews_N=',num2str(N),'_K=',num2str(K),'_M=',num2str(M),'_iid_',fn_property,'.mat'],'exp_set');

grey_rgb = (192/256)*ones(1,3);
dark_grey_rgb = (128/256)*ones(1,3);
figure
set(gca, 'FontSize', 14, 'LineWidth', 2)

xlim([0,200])
ylim([0,1.005])
xlabel('Solution ($x_0$)','interpreter','latex')
ylabel('Probability of $x_0$ in Subset','interpreter','latex')
%title(string_names{1},'interpreter','latex')

hold on
plot([1:200], (true_mean - min(true_mean))./(max(true_mean) - min(true_mean)), 'color', dark_grey_rgb, 'LineWidth', 1)
line([0,200], [1-alpha, 1-alpha], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5)
C1 = mean(S_testg_indicators_grad,2); %sum(all_S_indicators_d2(:,:,1),2)'/M;
C2 = mean(S_test2_indicators_grad,2); %sum(all_S_indicators_d2(:,:,1),2)'/M;
%scatter(feas_region(:,1), C, 'ko','markerfacecolor','k')
h3 = plot(feas_region(:,1), C4, 'r-', 'LineWidth', 2);
h4 = plot(feas_region(:,1), C1, 'g-', 'LineWidth', 2);
h5 = plot(feas_region(:,1), C2, 'b-', 'LineWidth', 2);
plot(exp_set,zeros(1,size(exp_set,2)),'kx','markerfacecolor','k', 'MarkerSize', 12, 'LineWidth', 1)

for i = 1:length(switch_on)
    x = [switch_on(i)-.5, switch_off(i)+.5, switch_off(i)+.5, switch_on(i)-.5];
    y = [0, 0, 1, 1]; 
    p=patch(x,y,'b','LineStyle','none','FaceAlpha',0.2);
end

legend([h3, h4, h5], 'GRPS', 'RPS w/ Dgrad w/ exact grad', 'RPS w/ D2 w/ exact grad')

box on
hold off


%%

figure

alpha = 0.05;
%load(['ctsnews_N=',num2str(N),'_K=',num2str(K),'_M=',num2str(M),'_iid_',fn_property,'.mat'],'exp_set');

grey_rgb = (192/256)*ones(1,3);
dark_grey_rgb = (128/256)*ones(1,3);
figure
set(gca, 'FontSize', 14, 'LineWidth', 2)

xlim([0,200])
ylim([0,1.005])
xlabel('Solution ($x_0$)','interpreter','latex')
ylabel('Probability of $x_0$ in Subset','interpreter','latex')
%title(string_names{1},'interpreter','latex')

hold on
plot([1:200], (true_mean - min(true_mean))./(max(true_mean) - min(true_mean)), 'color', dark_grey_rgb, 'LineWidth', 1)
line([0,200], [1-alpha, 1-alpha], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5)
%scatter(feas_region(:,1), C, 'ko','markerfacecolor','k')
h6 = plot(feas_region(:,1), C5, 'r-', 'LineWidth', 2);
h7 = plot(feas_region(:,1), C2, 'b-', 'LineWidth', 2);
plot(exp_set,zeros(1,size(exp_set,2)),'kx','markerfacecolor','k', 'MarkerSize', 12, 'LineWidth', 1)

for i = 1:length(switch_on)
    x = [switch_on(i)-.5, switch_off(i)+.5, switch_off(i)+.5, switch_on(i)-.5];
    y = [0, 0, 1, 1]; 
    p=patch(x,y,'b','LineStyle','none','FaceAlpha',0.2);
end

legend([h6, h7], 'RPS w/ d2', 'RPS w/ D2 w/ exact grad')

box on
hold off