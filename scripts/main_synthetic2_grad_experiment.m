% Script for testing the interface

clear;
clc;

add_rm_paths('add');

problem_string = 'synthetic2';
[oracle_string, oracle_n_rngs, feas_region, exp_set, k, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string] = init_problem(problem_string);

check_exceptions(discrep_string, fn_props, n_vec)

card_feas_region = size(feas_region, 1);
d = size(feas_region, 2);

opttol = 0.1; % delta

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
quiver(exp_set(:,1), exp_set(:,2), true_grad_exp_set(:,1), true_grad_exp_set(:,2), 'k:', 'LineWidth', 1.5) % plot true gradients

%title('True Performance Function and Gradients');
xlim([-2,2])
xlabel('$x^{(1)}$', 'Interpreter', 'latex')
ylim([-2,2])
ylabel('$x^{(2)}$', 'Interpreter', 'latex')
xticks([-2, -1, 0, 1, 2])
yticks([-2, -1, 0, 1, 2])
axis square
box on
hold off


%%
% Plot the ellipsoid of acceptable solutions
% delta=0.1 optimal = opttol

% (x1, y1) and (x2, y2) are coordinates of vertices on major axis
% e is eccentricity
% numbers from online ellipse calculator...
x1 = 1 - sqrt(10)/10;
x2 = 1 + sqrt(10)/10;
y1 = 1 - sqrt(10)/10;
y2 = 1 + sqrt(10)/10;
e = sqrt(6)/3;

a = 1/2*sqrt((x2-x1)^2+(y2-y1)^2);
b = a*sqrt(1-e^2);
t = linspace(0,2*pi);
X = a*cos(t);
Y = b*sin(t);
w = atan2(y2-y1,x2-x1);
xv = (x1+x2)/2 + X*cos(w) - Y*sin(w);
yv = (y1+y2)/2 + X*sin(w) + Y*cos(w);

hold on
h0 = plot(xv,yv,'k--', 'LineWidth', 2);
hold off

legend([h0], '$\mathcal{A}$', 'Location', 'northwest', 'Interpreter', 'latex')

% exportgraphics(gcf,'contours.png','Resolution',300)

%%

%%%% MAIN SECTION NOW (after first two subsections)

% if running multiple macroreplications...

% Number of macroreplications
%M = 1; % Small run to directly show subset
M = 100; % Larger run to show inclusion probabilities

% Calculate cutoffs for different methods
D_cutoff_PS = calc_cutoff(k, n_vec, alpha, 'ellinf'); % RPS: only performances
D_cutoff_PSG = calc_gradinf_cutoff(k, d, n_vec, alpha); % PSG
%D_cutoff_PSOG = calc_gradinf_cutoff(k, d-1, n_vec, alpha); % PSOG
D_cutoff_PSOG = calc_gradinftight_cutoff(k, d, n_vec, alpha); 

% Initialize data storage

% (Relaxed) Plausible Screening with dinf
S_PS_indicators = zeros(card_feas_region, M);
D_x0s = zeros(card_feas_region, M);
S_RPS_indicators = zeros(card_feas_region, M);

% Plausible Screening with Gradients
%S_indicators_gradinf = zeros(card_feas_region, M);
%D_x0s_gradinf = zeros(card_feas_region, M);
S_PSG_indicators = zeros(card_feas_region, M);

% Plausible Screening with Only Gradients
S_PSOG_indicators = zeros(card_feas_region, M);

print_problem_header(problem_string, feas_region, exp_set, fn_props)

parfor m = 1:M
    
    fprintf('Running macrorep %d of %d.\n', m, M)
    
    % SAMPLING
    
    % Generate data and calculate summary statistics
    fprintf('Generating sample data for plausible optima...\n')
    [sample_mean, sample_var, sample_mean_grad, sample_full_cov] = generate_grad_data(m, 'synthetic2_grad_oracle', oracle_n_rngs, exp_set, n_vec, 'grad');

%     % Original Plausible Screening
%     discrep_string = 'ellinf';
%     fprintf('Screening solutions for %s discrepancy...\n', discrep_string)
% %    [S_PS_indicators(:,m), D_x0s(:,m), S_RPS_indicators(:,m), ~] = PS_screen_fast(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ellinf', D_cutoff_PS, fn_props, prop_params, LP_solver_string, opttol);
%     S_RPS_indicators(:,m) = PS_screen_faster(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ellinf', D_cutoff_PS, fn_props, prop_params, LP_solver_string, opttol);
%     fprintf('\nResults of PS screening\n-------------------------------------------------------\n')
% %    fprintf('\t# of solutions in PS subset: \t\t\t%d\n', sum(S_PS_indicators(:,m)))  
%     fprintf('\t# of solutions in RPS subset: \t%d\n\n', sum(S_RPS_indicators(:,m)))

    % Gradient Plausible Screening w/ Dgradinf
    S_PSG_indicators(:,m) = PSG_screen_fast(feas_region, exp_set, sample_mean, sample_mean_grad, sample_full_cov, n_vec, D_cutoff_PSG, opttol);
    fprintf('\nResults of PSG screening (dinf) \n-------------------------------------------------------\n')
    fprintf('\t# of solutions in PSG subset: \t%d\n\n', sum(S_PSG_indicators(:,m)))
      
%     % Relaxed Gradient Plausible Screening w/ only gradients
%     S_PSOG_indicators(:,m) = PSOG_screen_fast(feas_region, exp_set, sample_mean, sample_mean_grad, sample_full_cov, n_vec, D_cutoff_PSOG, opttol);
%     fprintf('\nResults of PSOG screening \n-------------------------------------------------------\n')
%     fprintf('\t# of solutions in PSOG subset: \t%d\n\n', sum(S_PSOG_indicators(:,m)))

end

%% Setup for PS/PSG/PSOG limiting sets

[X_fine, Y_fine] = meshgrid(-1.975:.025:1.975, -1.975:.025:1.975);
n_grid = size(X_fine,1)*size(X_fine,2);
feas_region_fine = [reshape(X_fine,[n_grid, 1]), reshape(Y_fine,[n_grid, 1])];
card_feas_region_fine = size(feas_region_fine, 1);

x_vec_fine = feas_region_fine(:,1);
y_vec_fine = feas_region_fine(:,2);
true_mean_fine = x_vec_fine.^2 - x_vec_fine - x_vec_fine.*y_vec_fine - y_vec_fine + y_vec_fine.^2 + 1;
true_grad_fine = [2*x_vec_fine - y_vec_fine - 1, 2*y_vec_fine - x_vec_fine - 1];

x_mat_fine = reshape(x_vec_fine, [sqrt(card_feas_region_fine), sqrt(card_feas_region_fine)]);
y_mat_fine = reshape(y_vec_fine, [sqrt(card_feas_region_fine), sqrt(card_feas_region_fine)]);

%% Expensive setup for PS

%SOX_indicators = construct_det_no_grad_subset(feas_region, exp_set, true_mean, 'convex_nearopt', '', opttol);
SOX_indicators = construct_det_no_grad_subset(feas_region_fine, exp_set, true_mean_fine, 'convex_nearopt', '', opttol);


%% Plot test case for PS

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)
hold on

colormap summer
rescale_coef = max(max(true_mean_mat)) - min(min(true_mean_mat));
scatter(feas_region(:,1), feas_region(:,2), 20, mean(S_RPS_indicators,2)*rescale_coef, 's', 'filled')
contour(x_mat, y_mat, true_mean_mat, 10, 'LineColor', 'k')
scatter(exp_set(:,1), exp_set(:,2), 50, 'kx', 'LineWidth', 2) % plot experimental set X

h0 = plot(xv,yv,'k--', 'LineWidth', 2);
%SOX_indicators = construct_det_no_grad_subset(feas_region, exp_set, true_mean, 'convex_nearopt', '', opttol);
%[~,h1] = contour(x_mat, y_mat, reshape(SOX_indicators,[sqrt(card_feas_region), sqrt(card_feas_region)]), 1, 'LineColor', 'k', 'LineWidth', 2);
[~,h1] = contour(x_mat_fine, y_mat_fine, reshape(SOX_indicators,[sqrt(card_feas_region_fine), sqrt(card_feas_region_fine)]), 1, 'LineColor', 'k', 'LineWidth', 2);
legend([h1], '$\mathsf{S}^{\mathsf{O}}(\mathsf{X})$', 'Location', 'northwest', 'Interpreter', 'latex')

xlim([-2,2])
xlabel('$x^{(1)}$', 'Interpreter', 'latex')
ylim([-2,2])
ylabel('$x^{(2)}$', 'Interpreter', 'latex')
xticks([-2, -1, 0, 1, 2])
yticks([-2, -1, 0, 1, 2])
axis square
box on
hold off

% exportgraphics(gcf,'PS_heatmap.png','Resolution',300)

%% Plot test case for PSG

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)
hold on

colormap summer
rescale_coef = max(max(true_mean_mat)) - min(min(true_mean_mat));
scatter(feas_region(:,1), feas_region(:,2), 20, mean(S_PSG_indicators,2)*rescale_coef, 's', 'filled')
contour(x_mat, y_mat, true_mean_mat, 10, 'LineColor', 'k')
scatter(exp_set(:,1), exp_set(:,2), 50, 'kx', 'LineWidth', 2) % plot experimental set X

h0 = plot(xv, yv, 'k--', 'LineWidth', 2); % or patch
%SGX_indicators = construct_det_grad_subset(feas_region, exp_set, true_mean, true_grad, opttol);
%[~,h2] = contour(x_mat, y_mat, reshape(SGX_indicators,[sqrt(card_feas_region), sqrt(card_feas_region)]), 1, 'LineColor', 'k', 'LineWidth', 2);
SGX_indicators = construct_det_grad_subset(feas_region_fine, exp_set, true_mean_fine, true_grad_fine, opttol);
[~,h2] = contour(x_mat_fine, y_mat_fine, reshape(SGX_indicators,[sqrt(card_feas_region_fine), sqrt(card_feas_region_fine)]), 1, 'LineColor', 'k', 'LineWidth', 2);

legend([h2], '$\mathsf{S}^{\mathsf{G}}(\mathsf{X})$', 'Location', 'northwest', 'Interpreter', 'latex')


xlim([-2,2])
xlabel('$x^{(1)}$', 'Interpreter', 'latex')
ylim([-2,2])
ylabel('$x^{(2)}$', 'Interpreter', 'latex')
xticks([-2, -1, 0, 1, 2])
yticks([-2, -1, 0, 1, 2])
axis square
box on
hold off

% exportgraphics(gcf,'PSG_heatmap.png','Resolution',300)

%% Plot test case for PSOG

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)
hold on

colormap summer
rescale_coef = max(max(true_mean_mat)) - min(min(true_mean_mat));
scatter(feas_region(:,1), feas_region(:,2), 20, mean(S_PSOG_indicators,2)*rescale_coef, 's', 'filled')
contour(x_mat, y_mat, true_mean_mat, 10, 'LineColor', 'k')
scatter(exp_set(:,1), exp_set(:,2), 50, 'kx', 'LineWidth', 2) % plot experimental set X

h0 = plot(xv,yv,'k--', 'LineWidth', 2);
%SOGX_indicators = construct_det_grad_only_subset(feas_region, exp_set, true_mean, true_grad, opttol);
%[~,h3] = contour(x_mat, y_mat, reshape(SOGX_indicators,[sqrt(card_feas_region), sqrt(card_feas_region)]), 1, 'LineColor', 'k', 'LineWidth', 2);
SOGX_indicators = construct_det_grad_only_subset(feas_region_fine, exp_set, true_mean_fine, true_grad_fine, opttol);
[~,h3] = contour(x_mat_fine, y_mat_fine, reshape(SOGX_indicators,[sqrt(card_feas_region_fine), sqrt(card_feas_region_fine)]), 1, 'LineColor', 'k', 'LineWidth', 2);
legend([h3], '$\mathsf{S}^{\mathsf{OG}}(\mathsf{X})$', 'Location', 'northwest', 'Interpreter', 'latex')

xlim([-2,2])
xlabel('$x^{(1)}$', 'Interpreter', 'latex')
ylim([-2,2])
ylabel('$x^{(2)}$', 'Interpreter', 'latex')
xticks([-2, -1, 0, 1, 2])
yticks([-2, -1, 0, 1, 2])
axis square
box on
hold off

% exportgraphics(gcf,'PSOG_heatmap.png','Resolution',300)
