% Script for testing the interface

clear;
clc;

add_rm_paths('add');

problem_string = 'synthetic2';
[oracle_string, oracle_n_rngs, feas_region, exp_set, k, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string] = init_problem(problem_string);

check_exceptions(discrep_string, fn_props, n_vec)

card_feas_region = size(feas_region, 1);
d = size(feas_region, 2);

%%
% CALCULATE TRUE OBJECTIVE FUNCTIONS AND TRUE GRADIENTS

% mu(x) = x^2 - x - x*y - y + y^2 + 1
% GRAD mu(x) = [2x - y - 1, 2y - x - 1]
x_vec = feas_region(:,1);
y_vec = feas_region(:,2);
true_mean = x_vec.^2 - x_vec - x_vec.*y_vec - y_vec + y_vec.^2 + 1;
true_grad = [2*x_vec - y_vec - 1, 2*y_vec - x_vec - 1];
x_mat = reshape(x_vec, [sqrt(card_feas_region), sqrt(card_feas_region)]);
y_mat = reshape(y_vec, [sqrt(card_feas_region), sqrt(card_feas_region)]);
true_mean_mat = reshape(true_mean, [sqrt(card_feas_region), sqrt(card_feas_region)]);

% Extract values for exp_set
[~, exp_set_indices] = ismembertol(exp_set, feas_region, 'ByRows', true); 
true_mean_exp_set = true_mean(exp_set_indices,:);
true_grad_exp_set = true_grad(exp_set_indices,:);

%%
% PLOT TRUE FUNCTION/GRADIENTS/CONTOURS

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)
hold on

contour(x_mat, y_mat, true_mean_mat, 10, 'LineColor', 'k')
scatter(exp_set(:,1), exp_set(:,2), 50, 'kx', 'LineWidth', 2) % plot experimental set X
scatter(1, 1, 50, 'k*', 'LineWidth', 1) % plot optimal solution
quiver(exp_set(:,1), exp_set(:,2), true_grad_exp_set(:,1), true_grad_exp_set(:,2), 'b-', 'LineWidth', 1) % plot sample gradients

title('True Performance Function and Gradients');
xlim([-2,2])
xlabel('$x_1$', 'Interpreter', 'latex')
ylim([-2,2])
ylabel('$x_2$', 'Interpreter', 'latex')
axis square
box on
hold off


%%
% COMPUTE S(X) and plot it
SX_indicators = construct_det_subset(feas_region, exp_set, true_mean, 'convex', '');

% ...continuing from previous block of code
hold on
[~,h1] = contour(x_mat, y_mat, reshape(SX_indicators,[sqrt(card_feas_region), sqrt(card_feas_region)]), 1, 'LineColor', 'r', 'LineWidth', 2);
title('S(X) for PS with only performances')
hold off

%%
% COMPUTE SG(X) and plot it

SGX_indicators = construct_det_grad_subset(feas_region, exp_set, true_mean, true_grad);

% ...continuing from previous previous block of code
hold on
[~,h2] = contour(x_mat, y_mat, reshape(SGX_indicators,[sqrt(card_feas_region), sqrt(card_feas_region)]), 1, 'LineColor', 'm', 'LineWidth', 2);
title('SG(X) for GPS')
hold off

%%
% COMPUTE SRGO(X) and plot it
% GO stands for gradients only

SRGOX_indicators = construct_det_grad_only_subset(feas_region, exp_set, true_mean, true_grad);

% ...continuing from previous previous block of code
hold on
[~,h3] = contour(x_mat, y_mat, reshape(SRGOX_indicators,[sqrt(card_feas_region), sqrt(card_feas_region)]), 1, 'LineColor', 'c', 'LineWidth', 2);
title('SRGO(X) for GPS')
hold off

%%
% Add legends
legend([h1, h2, h3], 'S(X)', 'SG(X)', 'SRGO(X)', 'Location', 'southwest')
title('Deterministic Subsets')

%% RUN MACROREPLICATIONS

% Number of macroreplications
M = 1; % Small run to directly show subset
%M = 100; % Larger run to show inclusion probabilities

% Calculate cutoffs for PS
D_cutoff_d2 = calc_cutoff(k, n_vec, alpha, 'ell2'); % PS: only performances
D_cutoff_inf = calc_cutoff(k, n_vec, alpha, 'ellinf'); % PS: only performances
D_cutoff_grad = calc_grad_cutoff(k, d, n_vec, alpha); % GPS
D_cutoff_gradinf = calc_gradinf_cutoff(k, d, n_vec, alpha); % GPS
D_cutoff_gradinfeff = calc_gradinf_cutoff(k, d-1, n_vec, alpha); % only gradients

% Initialize data storage

% Plausible Screening
S_indicators_d2 = zeros(card_feas_region, M);
D_x0s = zeros(card_feas_region, M);
S_poly_indicators_d2 = zeros(card_feas_region, M);

% Gradient Plausible Screening with d2 grad and RGPS
S_indicators_grad = zeros(card_feas_region, M);
D_x0s_grad = zeros(card_feas_region, M);
S_poly_indicators_grad = zeros(card_feas_region, M);

% Gradient Plausible Screening with dinf grad and RGPS
S_indicators_gradinf = zeros(card_feas_region, M);
D_x0s_gradinf = zeros(card_feas_region, M);
S_poly_indicators_gradinf = zeros(card_feas_region, M);

% Relaxed relaxed PS, using performances too
Sonlygrad_d2_indicators = zeros(card_feas_region, M);
Sonlygrad_dinf_indicators = zeros(card_feas_region, M);

% Only gradients with dinf
Sonlygrad_dinf_eff_indicators = zeros(card_feas_region, M);

print_problem_header(problem_string, feas_region, exp_set, fn_props)

for m = 1:M
    
    fprintf('Running macrorep %d of %d.\n', m, M)
    
    % SAMPLING
    
    % Generate data and calculate summary statistics
    fprintf('Generating sample data for plausible optima...\n')
    [sample_mean, sample_var, sample_mean_grad, sample_full_cov] = generate_grad_data(m, 'synthetic2_grad_oracle', oracle_n_rngs, exp_set, n_vec, 'grad');

    % Original Plausible Screening
    discrep_string = 'ell2';
    fprintf('Screening solutions for %s discrepancy...\n', discrep_string)
    [S_indicators_d2(:,m), D_x0s(:,m), S_poly_indicators_d2(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ell2', D_cutoff_d2, fn_props, prop_params, LP_solver_string);
    fprintf('\nResults of PS screening\n---------------------------------------------feas----------\n')
    fprintf('\tstandardized discrepancy: \t\t\t\t%s\n', discrep_string)
    fprintf('\t# of solutions in PS d2 subset: \t\t\t%d\n', sum(S_indicators_d2(:,m)))  
    %fprintf('\t# of solutions in PS relaxed subset: \t%d\n\n', sum(S_poly_indicators_d2(:,m)))

    % Gradient Plausible Screening w/ Dgrad
    [S_indicators_grad(:,m), D_x0s_grad(:,m), S_poly_indicators_grad(:,m)] = GPS_screen(feas_region, exp_set, sample_mean, sample_mean_grad, sample_full_cov, n_vec, D_cutoff_grad);
    fprintf('\nResults of GPS screening\n-------------------------------------------------------\n')
    fprintf('\tstandardized discrepancy: \t\t\t\t Gradients\n')
    fprintf('\t# of solutions in GPS d2 subset: \t\t\t%d\n', sum(S_indicators_grad(:,m)))  
    %fprintf('\t# of solutions in GPS relaxed subset: \t%d\n\n', sum(S_poly_indicators_grad(:,m)))

    % Gradient Plausible Screening w/ Dgradinf
    [S_indicators_gradinf(:,m), D_x0s_gradinf(:,m), S_poly_indicators_gradinf(:,m)] = GPSinf_screen(feas_region, exp_set, sample_mean, sample_mean_grad, sample_full_cov, n_vec, D_cutoff_gradinf);
    fprintf('\nResults of GPS screening (dinf) \n-------------------------------------------------------\n')
    fprintf('\tstandardized discrepancy: \t\t\t\t Gradients\n')
    fprintf('\t# of solutions in GPS inf subset: \t\t\t%d\n', sum(S_indicators_gradinf(:,m)))  
    %fprintf('\t# of solutions in GPS inf relaxed subset: \t%d\n\n', sum(S_poly_indicators_gradinf(:,m)))
     
    % Relaxed Gradient Plausible Screening w/ only gradients
    [Sonlygrad_d2_indicators(:,m), Sonlygrad_dinf_indicators(:,m), Sonlygrad_dinf_eff_indicators(:,m)] = RGPS_onlygrad_screen(feas_region, exp_set, sample_mean, sample_mean_grad, sample_full_cov, n_vec, D_cutoff_grad, D_cutoff_gradinf, D_cutoff_gradinfeff);
    fprintf('\nResults of RGPS screening (only gradients) \n-------------------------------------------------------\n')
    %fprintf('\t# of solutions in RGPS (only gradients) d2 subset: \t\t\t%d\n', sum(Sonlygrad_d2_indicators(:,m)))  
    %fprintf('\t# of solutions in RGPS (only gradients) dinf subset: \t\t%d\n', sum(Sonlygrad_dinf_indicators(:,m)))
    fprintf('\t# of solutions in RGPS (only gradients) dinfeff subset: \t%d\n\n', sum(Sonlygrad_dinf_eff_indicators(:,m)))

end

%%
subsetstr = "onlygraddinfeff";

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)
hold on

% Print a subset from a single macroreplication
if strcmp(subsetstr, "nogradd2") == 1
    scatter(feas_region(S_indicators==1,1), feas_region(S_indicators==1,2), 'r.', 'filled')
elseif strcmp(subsetstr, "gradd2") == 1
    scatter(feas_region(S_indicators_grad==1,1), feas_region(S_indicators_grad==1,2), 'r.', 'filled')
    title('GPS with d2')
elseif strcmp(subsetstr, "graddinf") == 1 
    scatter(feas_region(S_indicators_gradinf==1,1), feas_region(S_indicators_gradinf==1,2), 'r.', 'filled')
    title('GPS with dinf')
elseif strcmp(subsetstr, "onlygraddinfeff") == 1
    scatter(feas_region(Sonlygrad_dinf_eff_indicators==1,1), feas_region(Sonlygrad_dinf_eff_indicators==1,2), 'r.', 'filled')
    title('RRGPS with dinf only gradients')
end

scatter(exp_set(:,1), exp_set(:,2), 50, 'kx', 'LineWidth', 2) % plot experimental set X
scatter(1, 1, 50, 'k*', 'LineWidth', 1) % plot optimal solution
quiver(exp_set(:,1), exp_set(:,2), sample_mean_grad(:,1), sample_mean_grad(:,2), 'b-', 'LineWidth', 1) % plot sample gradients
% divide MaxHeadSize by length of arrow to make arrow heads the same size
xlim([-2,2])
xlabel('$x_1$', 'Interpreter', 'latex')
ylim([-2,2])
ylabel('$x_2$', 'Interpreter', 'latex')
axis square
box on
hold off

%%
% if running multiple macroreplications...

% Number of macroreplications
%M = 1; % Small run to directly show subset
M = 100; % Larger run to show inclusion probabilities

% Calculate cutoffs for PS
D_cutoff_d2 = calc_cutoff(k, n_vec, alpha, 'ell2'); % PS: only performances
D_cutoff_inf = calc_cutoff(k, n_vec, alpha, 'ellinf'); % PS: only performances
D_cutoff_grad = calc_grad_cutoff(k, d, n_vec, alpha); % GPS
D_cutoff_gradinf = calc_gradinf_cutoff(k, d, n_vec, alpha); % GPS
D_cutoff_gradinfeff = calc_gradinf_cutoff(k, d-1, n_vec, alpha); % only gradients

% Initialize data storage

% Plausible Screening
S_indicators_d2 = zeros(card_feas_region, M);
D_x0s = zeros(card_feas_region, M);
S_poly_indicators_d2 = zeros(card_feas_region, M);

% Gradient Plausible Screening with d2 grad and RGPS
S_indicators_grad = zeros(card_feas_region, M);
D_x0s_grad = zeros(card_feas_region, M);
S_poly_indicators_grad = zeros(card_feas_region, M);

% Gradient Plausible Screening with dinf grad and RGPS
S_indicators_gradinf = zeros(card_feas_region, M);
D_x0s_gradinf = zeros(card_feas_region, M);
S_poly_indicators_gradinf = zeros(card_feas_region, M);

% Relaxed relaxed PS, using performances too
Sonlygrad_d2_indicators = zeros(card_feas_region, M);
Sonlygrad_dinf_indicators = zeros(card_feas_region, M);

% Only gradients with dinf
Sonlygrad_dinf_eff_indicators = zeros(card_feas_region, M);

print_problem_header(problem_string, feas_region, exp_set, fn_props)

parfor m = 1:M
    
    fprintf('Running macrorep %d of %d.\n', m, M)
    
    % SAMPLING
    
    % Generate data and calculate summary statistics
    fprintf('Generating sample data for plausible optima...\n')
    [sample_mean, sample_var, sample_mean_grad, sample_full_cov] = generate_grad_data(m, 'synthetic2_grad_oracle', oracle_n_rngs, exp_set, n_vec, 'grad');
% 
%     % Original Plausible Screening
%     discrep_string = 'ell2';
%     fprintf('Screening solutions for %s discrepancy...\n', discrep_string)
%     [S_indicators_d2(:,m), D_x0s(:,m), S_poly_indicators_d2(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ell2', D_cutoff_d2, fn_props, prop_params, LP_solver_string);
%     fprintf('\nResults of PS screening\n---------------------------------------------feas----------\n')
%     fprintf('\tstandardized discrepancy: \t\t\t\t%s\n', discrep_string)
%     fprintf('\t# of solutions in PS d2 subset: \t\t\t%d\n', sum(S_indicators_d2(:,m)))  
%     %fprintf('\t# of solutions in PS relaxed subset: \t%d\n\n', sum(S_poly_indicators_d2(:,m)))
% 
    % Gradient Plausible Screening w/ Dgrad
    [S_poly_indicators_grad(:,m)] = GPS_screen_fast(feas_region, exp_set, sample_mean, sample_mean_grad, sample_full_cov, n_vec, D_cutoff_grad);
    fprintf('\nResults of GPS screening\n-------------------------------------------------------\n')
    fprintf('\tstandardized discrepancy: \t\t\t\t Gradients\n')
    %fprintf('\t# of solutions in GPS d2 subset: \t\t\t%d\n', sum(S_indicators_grad(:,m)))  
    fprintf('\t# of solutions in GPS relaxed subset: \t%d\n\n', sum(S_poly_indicators_grad(:,m)))

%     % Gradient Plausible Screening w/ Dgradinf
%     [S_poly_indicators_gradinf(:,m)] = GPSinf_screen_fast(feas_region, exp_set, sample_mean, sample_mean_grad, sample_full_cov, n_vec, D_cutoff_gradinf);
%     fprintf('\nResults of GPS screening (dinf) \n-------------------------------------------------------\n')
%     fprintf('\tstandardized discrepancy: \t\t\t\t Gradients\n')
%     %fprintf('\t# of solutions in GPS inf subset: \t\t\t%d\n', sum(S_indicators_gradinf(:,m)))  
%     fprintf('\t# of solutions in GPS inf relaxed subset: \t%d\n\n', sum(S_poly_indicators_gradinf(:,m)))
%      
%     % Relaxed Gradient Plausible Screening w/ only gradients
%     [Sonlygrad_d2_indicators(:,m), Sonlygrad_dinf_indicators(:,m), Sonlygrad_dinf_eff_indicators(:,m)] = RGPS_onlygrad_screen(feas_region, exp_set, sample_mean, sample_mean_grad, sample_full_cov, n_vec, D_cutoff_grad, D_cutoff_gradinf, D_cutoff_gradinfeff);
%     fprintf('\nResults of RGPS screening (only gradients) \n-------------------------------------------------------\n')
%     %fprintf('\t# of solutions in RGPS (only gradients) d2 subset: \t\t\t%d\n', sum(Sonlygrad_d2_indicators(:,m)))  
%     %fprintf('\t# of solutions in RGPS (only gradients) dinf subset: \t\t%d\n', sum(Sonlygrad_dinf_indicators(:,m)))
%     fprintf('\t# of solutions in RGPS (only gradients) dinfeff subset: \t%d\n\n', sum(Sonlygrad_dinf_eff_indicators(:,m)))

end

%% Plot test case

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)
hold on

colormap summer
rescale_coef = max(max(true_mean_mat)) - min(min(true_mean_mat));
scatter(feas_region(:,1), feas_region(:,2), 20, mean(Sonlygrad_dinf_eff_indicators,2)*rescale_coef, 's', 'filled')
%scatter(feas_region(:,1), feas_region(:,2), 20, mean(Sonlygrad_dinf_eff_indicators,2), 's', 'filled')
contour(x_mat, y_mat, true_mean_mat, 10, 'LineColor', 'k')
scatter(exp_set(:,1), exp_set(:,2), 50, 'kx', 'LineWidth', 2) % plot experimental set X
scatter(1, 1, 50, 'k*', 'LineWidth', 1) % plot optimal solution

title('Inclusion Probabilities for RRGPS');
xlim([-2,2])
xlabel('$x_1$', 'Interpreter', 'latex')
ylim([-2,2])
ylabel('$x_2$', 'Interpreter', 'latex')
axis square
box on
hold off

%% Plot test case

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)
hold on

colormap summer
rescale_coef = max(max(true_mean_mat)) - min(min(true_mean_mat));
scatter(feas_region(:,1), feas_region(:,2), 20, mean(S_poly_indicators_gradinf,2)*rescale_coef, 's', 'filled')
%scatter(feas_region(:,1), feas_region(:,2), 20, mean(Sonlygrad_dinf_eff_indicators,2), 's', 'filled')
contour(x_mat, y_mat, true_mean_mat, 10, 'LineColor', 'k')
scatter(exp_set(:,1), exp_set(:,2), 50, 'kx', 'LineWidth', 2) % plot experimental set X
scatter(1, 1, 50, 'k*', 'LineWidth', 1) % plot optimal solution

title('Inclusion Probabilities for RGPS dinf');
xlim([-2,2])
xlabel('$x_1$', 'Interpreter', 'latex')
ylim([-2,2])
ylabel('$x_2$', 'Interpreter', 'latex')
axis square
box on
hold off


%% Plot test case

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)
hold on

colormap summer
rescale_coef = max(max(true_mean_mat)) - min(min(true_mean_mat));
scatter(feas_region(:,1), feas_region(:,2), 20, mean(S_poly_indicators_grad,2)*rescale_coef, 's', 'filled')
%scatter(feas_region(:,1), feas_region(:,2), 20, mean(Sonlygrad_dinf_eff_indicators,2), 's', 'filled')
contour(x_mat, y_mat, true_mean_mat, 10, 'LineColor', 'k')
scatter(exp_set(:,1), exp_set(:,2), 50, 'kx', 'LineWidth', 2) % plot experimental set X
scatter(1, 1, 50, 'k*', 'LineWidth', 1) % plot optimal solution

title('Inclusion Probabilities for RGPS d2');
xlim([-2,2])
xlabel('$x_1$', 'Interpreter', 'latex')
ylim([-2,2])
ylabel('$x_2$', 'Interpreter', 'latex')
axis square
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
