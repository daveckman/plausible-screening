% Make plots for tandem experiment

%% Sorted discrepancy plot

% Load shared data
myVars = {'D_cutoffs','card_feas_region'};
load('tandem_M=1_iid_ellinf_convex_budget50.mat', myVars{:});

% Load standardized discrepancies
load('tandem_M=1_iid_ell1_convex_budget50.mat', 'D_x0s');
D_x0s_d1 = D_x0s;
load('tandem_M=1_iid_ell2_convex_budget50.mat', 'D_x0s');
D_x0s_d2 = D_x0s;
load('tandem_M=1_iid_ellinf_convex_budget50.mat', 'D_x0s');
D_x0s_dinf = D_x0s;


% Plotting Setup
string_names = {'PO: $d^1$', 'PO: $d^2$', 'PO: $d^{\infty}$'};
rgb_red = [0.8500, 0.3250, 0.0980]; 
rgb_yellow = [0.9290, 0.6940, 0.1250];
rgb_purple = [0.4940, 0.1840, 0.5560];

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)
xlim([0, card_feas_region])
xlabel('Sorted Solutions','interpreter','latex')
ylabel('$\log(D(x_0, \widehat{\mu}, \widehat{\Sigma}, n)/\mathsf{D})$','interpreter','latex')
hold on
%plot(1:card_feas_region, sort_D_x0s, 'b-', 'LineWidth', 2);
h1 = plot(log(sort(D_x0s_d1)/D_cutoffs(1)), 'color', rgb_red, 'LineWidth', 2);
h2 = plot(log(sort(D_x0s_d2)/D_cutoffs(2)), 'color', rgb_yellow, 'LineWidth', 2);
h3 = plot(log(sort(D_x0s_dinf)/D_cutoffs(3)), 'color', rgb_purple, 'LineWidth', 2);
line([0,card_feas_region], [0, 0], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5);
line([sum(D_x0s_d1 <= D_cutoffs(1)), sum(D_x0s_d1 <= D_cutoffs(1))], [-2, 0], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5);
line([sum(D_x0s_d2 <= D_cutoffs(2)), sum(D_x0s_d2 <= D_cutoffs(2))], [-2, 0], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5);
line([sum(D_x0s_dinf <= D_cutoffs(3)), sum(D_x0s_dinf <= D_cutoffs(3))], [-2, 0], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5);
%plot([0, card_feas_region], [0, 0], 'k:', 'LineWidth', 1.5);
legend([h1, h2, h3], string_names, 'location', 'northwest', 'Interpreter', 'latex');
legend boxoff
txt = 'Cutoff';
text(225000, 0.25, txt, 'FontSize', 14)
box on
hold off

%print(['sorted_min_discrep_tandem_all'],'-dpng','-r500')

%% Stacked histogram

load('tandem_true_mean_data_n=500_crn.mat', 'sample_mean')
load('tandem_M=1_iid_ell2_convex_budget50.mat', 'S_indicators')
load('tandem_M=1_iid_ell2_mode=2_convex_amended.mat', 'S_poly_indicators')

opt_gap = sample_mean - min(sample_mean);
max_opt_gap = max(opt_gap);
%bins = linspace(0, max_opt_gap, 20);
bins = 0:5:100;

% Calculate histogram data
[raw_X,~] = histcounts(opt_gap, bins);
[raw_Spoly,~] = histcounts(opt_gap(S_poly_indicators == 1), bins);
[raw_S,~] = histcounts(opt_gap(S_indicators == 1), bins);
bin_centers = 2.5:5:97.5;

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)
xlim([0, 100])
xlabel('Optimality Gap','interpreter','latex')
ylabel('Number of Solutions','interpreter','latex')

hold on
bar(bin_centers, [raw_S; raw_Spoly - raw_S; raw_X - raw_Spoly]', 1, 'stacked', 'LineWidth', 1)
%histogram(opt_gap, bins);
%histogram(opt_gap(S_poly_indicators == 1), bins);
%histogram(opt_gap(S_indicators == 1), bins);

legend_strings = {'$\mathcal{S}$', '$\mathcal{S}_{\mathrm{poly}}$', '$\mathcal{X}$'};
legend(legend_strings, 'location', 'northeast', 'Interpreter', 'latex');
legend boxoff
hold off

print('histogram_opt_gap_tandem','-dpng','-r500')