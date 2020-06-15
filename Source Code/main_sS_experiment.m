% Script for testing the interface

clear;
clc;

% SETUP
% Run a script to initialize:
%   oracle_string: string for the name of the oracle function
%   oracle_n_rngs: number of rngs needed by the oracle
%   feas_region: feasible region matrix (each row corresponds to a solution)
%   exp_set: experimental set matrix (each row corresponds to a solution)
%   k: number of evaluated solutions
%   n_vec: column vector of sample sizes for each evaluated solution
%   alpha: confidence level = 1-alpha
%   discrep_string: string for discrepancy type {'ell1', 'ell2', 'ellinf', 'CRN'}
%   fn_props: string for functional property {'convex', 'lipschitz', 'lipschitz_proj}
%   prop_params: Lipschitz constant (if applicable)

test_sS_inventory;


%%
% CHECK FOR EXCEPTIONS

accept_discreps = {'ell1', 'ell2', 'ellinf', 'CRN'};
if ~any(strcmp(discrep_string, accept_discreps))
    fprintf('\nERROR: "%s" is not a valid standardized discrepancy.\n', discrep_string)
    fprintf('Please specify a valid standardized discrepancy:')
    fprintf('\t%s', accept_discreps{:})
    fprintf('.\n')
    return
end

accept_fn_props = {'lipschitz', 'lipschitz_proj', 'convex'};
if ~any(strcmp(fn_props, accept_fn_props))
    fprintf('\nERROR: "%s" is not a valid functional property.\n', fn_props)
    fprintf('Please specify a valid functional property:')
    fprintf('\t%s', accept_fn_props{:})
    fprintf('.\n')
    return
end

if strcmp(discrep_string, 'CRN') == 1
    % Check that all values in n_vec vector are equal
    if min(n_vec) ~= max(n_vec)
        fprintf('All sample sizes must be equal when using CRN.\n')
        return
    end
    
    % Check if n >= k so that sample_var is non-singular
    if n_vec(1) < k
        fprintf('Common sample size n = %d must be at least k = %d.\n', n_vec(1), k)
        return
    end
end

%% RUN MACROREPLICATIONS

M = 100; % Number of macroreplications
card_feas_region = size(feas_region, 1);

% Initialize data storage
S_indicators_d1 = zeros(card_feas_region, M);
S_indicators_d2 = zeros(card_feas_region, M);
S_indicators_dinf = zeros(card_feas_region, M);
%S_indicators_dcrn = zeros(card_feas_region, M);
S_poly_indicators_d1 = zeros(card_feas_region, M);
S_poly_indicators_d2 = zeros(card_feas_region, M);
S_poly_indicators_dinf = zeros(card_feas_region, M);
%S_poly_indicators_dcrn = zeros(card_feas_region, M);


for m = 1:M
    
    fprintf('Running macrorep %d of %d.\n', m, M)
    
    % SAMPLING
    
    % Generate data and calculate summary statistics
    fprintf('Generating sample data in parallel...\n')
    [sample_mean, sample_var] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, discrep_string);

    % SCREENING (USING DIFFERENT DISCREPANCIES)
    discrep_string = 'ell1';
    fprintf('Screening solutions in parallel for %s discrepancy...\n', discrep_string)
    [S_indicators_d1(:,m), ~, S_poly_indicators_d1(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string);
    fprintf('\nResults of PO screening\n-------------------------------------------------------\n')
    fprintf('\tstandardized discrepancy: \t\t\t\t%s\n', discrep_string)
    fprintf('\t# of solutions in PO subset: \t\t\t%d\n', sum(S_indicators_d1(:,m)))
    fprintf('\t# of solutions in PO relaxed subset: \t%d\n\n', sum(S_poly_indicators_d1(:,m)))

    discrep_string = 'ell2';
    fprintf('Screening solutions in parallel for %s discrepancy...\n', discrep_string)
    [S_indicators_d2(:,m), ~, S_poly_indicators_d2(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string);
    fprintf('\nResults of PO screening\n-------------------------------------------------------\n')
    fprintf('\tstandardized discrepancy: \t\t\t\t%s\n', discrep_string)
    fprintf('\t# of solutions in PO subset: \t\t\t%d\n', sum(S_indicators_d2(:,m)))  
    fprintf('\t# of solutions in PO relaxed subset: \t%d\n\n', sum(S_poly_indicators_d2(:,m)))

    discrep_string = 'ellinf';
    fprintf('Screening solutions in parallel for %s discrepancy...\n', discrep_string)
    [S_indicators_dinf(:,m), ~, S_poly_indicators_dinf(:,m), ~] = PO_screen(feas_region, exp_set, sample_mean, sample_var, n_vec, alpha, discrep_string, fn_props, prop_params, LP_solver_string);
    fprintf('\nResults of PO screening\n-------------------------------------------------------\n')
    fprintf('\tstandardized discrepancy: \t\t\t\t%s\n', discrep_string)
    fprintf('\t# of solutions in PO subset: \t\t\t%d\n', sum(S_indicators_dinf(:,m)))
    fprintf('\t# of solutions in PO relaxed subset: \t%d\n\n', sum(S_poly_indicators_dinf(:,m)))

end

%%
% PLOTTING SUBSETS

figure
subplot(1, 2, 1);

C = sum(S_indicators_d1,2)'/M;
for k = 1:size(feas_region,1)
    rectangle('Position',[feas_region(k,1)-C(k)/2  feas_region(k,2)-C(k)/2 C(k) C(k)],'facecolor',[0.2 0.8 0.8],'linestyle','none','Curvature',1)
end
axis square
hold on

plot(exp_set(:,1),exp_set(:,2),'kd','markerfacecolor','k')
vch = convhull(feas_region(:,1),feas_region(:,2));
plot(feas_region(vch,1),feas_region(vch,2),'k-')
xlim([min(feas_region(:,1))-5,max(feas_region(:,1))+5])
ylim([min(feas_region(:,2))-5,max(feas_region(:,2))+5])
xlabel('Reorder quantity','interpreter','latex')
ylabel('Order-to quantity','interpreter','latex')
title('$d^1$ discrepancy (standard PO)', 'interpreter', 'latex')
hold off

subplot(1, 2, 2);

C = sum(S_poly_indicators_d1,2)'/M;
for k = 1:size(feas_region,1)
    rectangle('Position',[feas_region(k,1)-C(k)/2  feas_region(k,2)-C(k)/2 C(k) C(k)],'facecolor',[0.8 0.2 0.2],'linestyle','none','Curvature',1)
end
axis square
hold on

plot(exp_set(:,1),exp_set(:,2),'kd','markerfacecolor','k')
vch = convhull(feas_region(:,1),feas_region(:,2));
plot(feas_region(vch,1),feas_region(vch,2),'k-')
xlim([min(feas_region(:,1))-5,max(feas_region(:,1))+5])
ylim([min(feas_region(:,2))-5,max(feas_region(:,2))+5])
xlabel('Reorder quantity','interpreter','latex')
ylabel('Order-to quantity','interpreter','latex')
title('$d^1$ discrepancy (relaxed PO)', 'interpreter', 'latex')
hold off

figure;
subplot(1, 2, 1);

C = sum(S_indicators_d2,2)'/M;
for k = 1:size(feas_region,1)
    rectangle('Position',[feas_region(k,1)-C(k)/2  feas_region(k,2)-C(k)/2 C(k) C(k)],'facecolor',[0.2 0.8 0.8],'linestyle','none','Curvature',1)
end
axis square
hold on

plot(exp_set(:,1),exp_set(:,2),'kd','markerfacecolor','k')
vch = convhull(feas_region(:,1),feas_region(:,2));
plot(feas_region(vch,1),feas_region(vch,2),'k-')
xlim([min(feas_region(:,1))-5,max(feas_region(:,1))+5])
ylim([min(feas_region(:,2))-5,max(feas_region(:,2))+5])
xlabel('Reorder quantity','interpreter','latex')
ylabel('Order-to quantity','interpreter','latex')
title('$d^2$ discrepancy (standard PO)', 'interpreter', 'latex')
hold off

subplot(1, 2, 2);

C = sum(S_poly_indicators_d2,2)'/M;
for k = 1:size(feas_region,1)
    rectangle('Position',[feas_region(k,1)-C(k)/2  feas_region(k,2)-C(k)/2 C(k) C(k)],'facecolor',[0.8 0.2 0.2],'linestyle','none','Curvature',1)
end
axis square
hold on

plot(exp_set(:,1),exp_set(:,2),'kd','markerfacecolor','k')
vch = convhull(feas_region(:,1),feas_region(:,2));
plot(feas_region(vch,1),feas_region(vch,2),'k-')
xlim([min(feas_region(:,1))-5,max(feas_region(:,1))+5])
ylim([min(feas_region(:,2))-5,max(feas_region(:,2))+5])
xlabel('Reorder quantity','interpreter','latex')
ylabel('Order-to quantity','interpreter','latex')
title('$d^2$ discrepancy (relaxed PO)', 'interpreter', 'latex')
hold off

figure;
subplot(1, 2, 1)

C = sum(S_indicators_dinf,2)'/M;
for k = 1:size(feas_region,1)
    rectangle('Position',[feas_region(k,1)-C(k)/2  feas_region(k,2)-C(k)/2 C(k) C(k)],'facecolor',[0.2 0.8 0.8],'linestyle','none','Curvature',1)
end
axis square
hold on

plot(exp_set(:,1),exp_set(:,2),'kd','markerfacecolor','k')
vch = convhull(feas_region(:,1),feas_region(:,2));
plot(feas_region(vch,1),feas_region(vch,2),'k-')
xlim([min(feas_region(:,1))-5,max(feas_region(:,1))+5])
ylim([min(feas_region(:,2))-5,max(feas_region(:,2))+5])
xlabel('Reorder quantity','interpreter','latex')
ylabel('Order-to quantity','interpreter','latex')
title('$d^\infty$ discrepancy (standard PO)', 'interpreter', 'latex')
hold off

subplot(1, 2, 2);

C = sum(S_poly_indicators_dinf,2)'/M;
for k = 1:size(feas_region,1)
    rectangle('Position',[feas_region(k,1)-C(k)/2  feas_region(k,2)-C(k)/2 C(k) C(k)],'facecolor',[0.8 0.2 0.2],'linestyle','none','Curvature',1)
end
axis square
hold on

plot(exp_set(:,1),exp_set(:,2),'kd','markerfacecolor','k')
vch = convhull(feas_region(:,1),feas_region(:,2));
plot(feas_region(vch,1),feas_region(vch,2),'k-')
xlim([min(feas_region(:,1))-5,max(feas_region(:,1))+5])
ylim([min(feas_region(:,2))-5,max(feas_region(:,2))+5])
xlabel('Reorder quantity','interpreter','latex')
ylabel('Order-to quantity','interpreter','latex')
title('$d^\infty$ discrepancy (relaxed PO)', 'interpreter', 'latex')
hold off

%% PLOTTING SUBSET SIZES (ECDFS)

% Compute subset sizes
subsize_d1 = sum(S_indicators_d1, 1);
subsize_d1_poly = sum(S_poly_indicators_d1, 1);
subsize_d2 = sum(S_indicators_d2, 1);
subsize_d2_poly = sum(S_poly_indicators_d2, 1);
subsize_dinf = sum(S_indicators_dinf, 1);
subsize_dinf_poly = sum(S_poly_indicators_dinf, 1);

% Compute empirical cdf of expected sample size
ecdf_d1 = zeros(1, card_feas_region);
ecdf_d1_poly = zeros(1, card_feas_region);
ecdf_d2 = zeros(1, card_feas_region);
ecdf_d2_poly = zeros(1, card_feas_region);
ecdf_dinf = zeros(1, card_feas_region);
ecdf_dinf_poly = zeros(1, card_feas_region);

for j = 1:card_feas_region
    ecdf_d1(j) = mean(subsize_d1 <= j);
    ecdf_d1_poly(j) = mean(subsize_d1_poly <= j);
    ecdf_d2(j) = mean(subsize_d2 <= j);
    ecdf_d2_poly(j) = mean(subsize_d2_poly <= j);
    ecdf_dinf(j) = mean(subsize_dinf <= j);
    ecdf_dinf_poly(j) = mean(subsize_dinf_poly <= j);
end

figure
hold on

plot([1:card_feas_region], ecdf_d1, 'b-', 'LineWidth', 2);
plot([1:card_feas_region], ecdf_d1_poly, 'b:','LineWidth', 2);
plot([1:card_feas_region], ecdf_d2, 'r-', 'LineWidth', 2);
plot([1:card_feas_region], ecdf_d2_poly, 'r:', 'LineWidth', 2);
plot([1:card_feas_region], ecdf_dinf, 'g-', 'LineWidth', 2);
plot([1:card_feas_region], ecdf_dinf_poly, 'g:', 'LineWidth', 2);

xlabel('Subset Size ($s$)', 'Interpreter', 'latex')
ylabel('P($|S| \leq s$)', 'Interpreter', 'latex')
leg1 = legend('$d^1$ standard', '$d^1$ relaxed', '$d^2$ standard', '$d^2$ relaxed', '$d^{\infty}$ standard', '$d^{\infty}$ relaxed', 'location', 'southeast');
legend boxoff
set(leg1,'Interpreter','latex');
set(gca, 'FontSize', 14, 'LineWidth', 2)
set(gcf,'units','inches','position',[1,1,9,4.5])

hold off