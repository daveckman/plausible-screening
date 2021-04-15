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
quiver(exp_set(:,1), exp_set(:,2), true_grad_exp_set(:,1), true_grad_exp_set(:,2), 'b-', 'LineWidth', 1) % plot sample gradients

%title('True Performance Function and Gradients');
xlim([-2,2])
xlabel('$x^{(1)}$', 'Interpreter', 'latex')
ylim([-2,2])
ylabel('$x^{(2)}$', 'Interpreter', 'latex')
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
h0 = plot(xv,yv,'b-', 'LineWidth', 2);
hold off

%%
% COMPUTE SO(X), SG(X), and SOG(X) and plot them
SOX_indicators = construct_det_no_grad_subset(feas_region, exp_set, true_mean, 'convex_nearopt', '', opttol);
SGX_indicators = construct_det_grad_subset(feas_region, exp_set, true_mean, true_grad, opttol);
SOGX_indicators = construct_det_grad_only_subset(feas_region, exp_set, true_mean, true_grad, opttol);

% ...continuing from previous block of code
hold on
[~,h1] = contour(x_mat, y_mat, reshape(SOX_indicators,[sqrt(card_feas_region), sqrt(card_feas_region)]), 1, 'LineColor', 'r', 'LineWidth', 2);
[~,h2] = contour(x_mat, y_mat, reshape(SGX_indicators,[sqrt(card_feas_region), sqrt(card_feas_region)]), 1, 'LineColor', 'm', 'LineWidth', 2);
[~,h3] = contour(x_mat, y_mat, reshape(SOGX_indicators,[sqrt(card_feas_region), sqrt(card_feas_region)]), 1, 'LineColor', 'c', 'LineWidth', 2);
hold off

legend([h0, h1, h2, h3], '$\mathcal{A}$', '$\mathrm{S}^{\mathrm{O}}(\mathrm{X})$', '$\mathrm{S}^{\mathrm{G}}(\mathrm{X})$', '$\mathrm{S}^{\mathrm{OG}}(\mathrm{X})$', 'Location', 'northwest', 'Interpreter', 'latex')

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

%%%% MAIN SECTION NOW (after first two subsections)

% if running multiple macroreplications...

% Number of macroreplications
M = 1; % Small run to directly show subset
%M = 100; % Larger run to show inclusion probabilities

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

    % Original Plausible Screening
    discrep_string = 'ellinf';
    fprintf('Screening solutions for %s discrepancy...\n', discrep_string)
    [S_PS_indicators(:,m), D_x0s(:,m), S_RPS_indicators(:,m), ~] = PS_screen_fast(feas_region, exp_set, sample_mean, sample_var, n_vec, 'ellinf', D_cutoff_PS, fn_props, prop_params, LP_solver_string, opttol);
    fprintf('\nResults of PS screening\n-------------------------------------------------------\n')
    %fprintf('\tstandardized discrepancy: \t\t\t\t%s\n', discrep_string)
    fprintf('\t# of solutions in PS subset: \t\t\t%d\n', sum(S_PS_indicators(:,m)))  
    fprintf('\t# of solutions in RPS subset: \t%d\n\n', sum(S_RPS_indicators(:,m)))


%     % Gradient Plausible Screening w/ Dgradinf
%     S_PSG_indicators(:,m) = PSG_screen_fast(feas_region, exp_set, sample_mean, sample_mean_grad, sample_full_cov, n_vec, D_cutoff_PSG, opttol);
%     fprintf('\nResults of PSG screening (dinf) \n-------------------------------------------------------\n')
%     %fprintf('\tstandardized discrepancy: \t\t\t\t Gradients\n')
%     fprintf('\t# of solutions in PSG subset: \t%d\n\n', sum(S_PSG_indicators(:,m)))
      
%     % Relaxed Gradient Plausible Screening w/ only gradients
%     S_PSOG_indicators(:,m) = PSOG_screen_fast(feas_region, exp_set, sample_mean, sample_mean_grad, sample_full_cov, n_vec, D_cutoff_PSOG, opttol);
%     fprintf('\nResults of PSOG screening \n-------------------------------------------------------\n')
%     %fprintf('\tstandardized discrepancy: \t\t\t\t Gradients\n')
%     fprintf('\t# of solutions in PSOG subset: \t%d\n\n', sum(S_PSOG_indicators(:,m)))

end

%% Plot test case

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)
hold on

colormap summer
rescale_coef = max(max(true_mean_mat)) - min(min(true_mean_mat));
scatter(feas_region(:,1), feas_region(:,2), 20, mean(S_RPS_indicators,2)*rescale_coef, 's', 'filled')
contour(x_mat, y_mat, true_mean_mat, 10, 'LineColor', 'k')
scatter(exp_set(:,1), exp_set(:,2), 50, 'kx', 'LineWidth', 2) % plot experimental set X
scatter(1, 1, 50, 'k*', 'LineWidth', 1) % plot optimal solution

%title('Inclusion Probabilities for PS');
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
scatter(feas_region(:,1), feas_region(:,2), 20, mean(S_PSG_indicators,2)*rescale_coef, 's', 'filled')
contour(x_mat, y_mat, true_mean_mat, 10, 'LineColor', 'k')
scatter(exp_set(:,1), exp_set(:,2), 50, 'kx', 'LineWidth', 2) % plot experimental set X
scatter(1, 1, 50, 'k*', 'LineWidth', 1) % plot optimal solution

%title('Inclusion Probabilities for PSG');
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
scatter(feas_region(:,1), feas_region(:,2), 20, mean(S_PSOG_indicators,2)*rescale_coef, 's', 'filled')
contour(x_mat, y_mat, true_mean_mat, 10, 'LineColor', 'k')
scatter(exp_set(:,1), exp_set(:,2), 50, 'kx', 'LineWidth', 2) % plot experimental set X
scatter(1, 1, 50, 'k*', 'LineWidth', 1) % plot optimal solution

%title('Inclusion Probabilities for PSOG');
xlim([-2,2])
xlabel('$x_1$', 'Interpreter', 'latex')
ylim([-2,2])
ylabel('$x_2$', 'Interpreter', 'latex')
axis square
box on
hold off


%% THE REST IS NOT USED

% %% 
% % CHECK RGPS WITH PLUGGING IN TRUE GRADIENTS
% %M = 100;
% 
% exact_grads = true_grad(exp_set,:);
% 
% D_grad = calc_grad_cutoff(k, d, n_vec, alpha);
% D_d2 = calc_cutoff(k, n_vec, alpha, 'ell2');
% 
% S_test_indicators_grad = zeros(card_feas_region, M);
% 
% parfor m = 1:M
%     
%     fprintf('Running macrorep %d of %d.\n', m, M)
%     
%     % SAMPLING
%     
%     % Generate data and calculate summary statistics
%     fprintf('Generating sample data for plausible optima...\n')
%     [sample_mean, sample_var, sample_mean_grad, sample_full_cov] = generate_grad_data(m, 'cts_newsvendor_grad_oracle', oracle_n_rngs, exp_set, n_vec, 'grad');
% 
%     
%     
%     %[sample_mean, sample_var] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, 'ell1');
%     S_testg_indicators_grad(:,m) = RGPSexact_screen(feas_region, exp_set, sample_mean, sample_var, exact_grads, n_vec, D_grad);
%     S_test2_indicators_grad(:,m) = RGPSexact_screen(feas_region, exp_set, sample_mean, sample_var, exact_grads, n_vec, D_d2);
%     fprintf('\nResults of RGPS screening (using exact gradients)\n-------------------------------------------------------\n')
%     fprintf('\tstandardized discrepancy: \t\t\t\t Gradients\n')
%     fprintf('\t# of solutions in GPS relaxed subset: \t%d\n\n', sum(S_testg_indicators_grad(:,m)))
%     fprintf('\t# of solutions in GPS relaxed subset (d2 cutoff): \t%d\n\n', sum(S_test2_indicators_grad(:,m)))
% 
% end
