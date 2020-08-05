clear
clc;

fn_property = 'lipschitz_proj';
%fn_property = 'lipschitz_proj'; % 'convex' or 'lipschitz_proj'

if strcmp(fn_property, 'lipschitz_proj') == 1
    plt_title = 'Lipschitz';
elseif strcmp(fn_property, 'convex') == 1
    plt_title = 'Convex';
end

setup = 'fixedN'; % 'fixedK' or 'fixedN'

if strcmp(setup, 'fixedN') == 1
    N_list = [400, 400, 400, 400];
    K_list = [5, 10, 20, 50];
    x_axis_limits = [min(K_list), max(K_list)];
    x_axis_pts = K_list;
    x_axis_label = 'Cardinality of Experimental Set ($k$)';
elseif strcmp(setup, 'fixedK') == 1
    N_list = [400, 600, 1000, 2000, 4000];
    K_list = [5, 5, 5, 5, 5];
    x_axis_limits = [min(N_list), max(N_list)];
    x_axis_pts = N_list;
    x_axis_label = 'Total Sample Size';
end

%% Test plotting for each experiment
% 
% K_list = [5, 10, 20, 50];
% 
% myVars = {'SS_indicators', 'S_indicators_d1', 'S_poly_indicators_d1', 'S_indicators_d2', 'S_poly_indicators_d2', 'S_indicators_dinf', 'S_poly_indicators_dinf'};
% string_names = {'ExtSTB', 'PO: $d^1$', 'PO (relaxed): $d^1$', 'PO: $d^2$', 'PO (relaxed): $d^2$', 'PO: $d^{\infty}$', 'PO (relaxed): $d^{\infty}$'};
% colors = {'k:', 'b-', 'b:', 'g-', 'g:', 'm-', 'm:'};
%     
% for i = 1:length(K_list)
%     load(['ctsnews_N=400_K=',num2str(K_list(i)),'_iid_lipschitz.mat'],myVars{:});
%     all_S_indicators = {SS_indicators, S_indicators_d1, S_poly_indicators_d1, S_indicators_d2, S_poly_indicators_d2, S_indicators_dinf, S_poly_indicators_dinf};
%     
%     plot_sample_size_ecdfs(all_S_indicators, string_names, colors)
% end

%% Experiments for fixed N (total number of replications)
myVars = {'SS_indicators', 'S_indicators_d1', 'S_indicators_d2', 'S_indicators_dinf'};

% Initialize
card_feas_region = 200;
macroreps = 3000;
all_SS_indicators = zeros(card_feas_region, macroreps, length(K_list));
all_S_indicators_d1 = zeros(card_feas_region, macroreps, length(K_list));
all_S_indicators_d2 = zeros(card_feas_region, macroreps, length(K_list));
all_S_indicators_dinf = zeros(card_feas_region, macroreps, length(K_list));

% Read data into matrix
for i = 1:length(K_list)
    load(['ctsnews_N=',num2str(N_list(i)),'_K=',num2str(K_list(i)),'_M=3000_iid_',fn_property,'.mat'],myVars{:});
    all_SS_indicators(:,:,i) = SS_indicators;
    all_S_indicators_d1(:,:,i) = S_indicators_d1;
    all_S_indicators_d2(:,:,i) = S_indicators_d2;
    all_S_indicators_dinf(:,:,i) = S_indicators_dinf;
end

%% Compute summary statistics

% Compute P(x0 in S)
inc_probs_SS = mean(all_SS_indicators,2);
inc_probs_S_d1 = mean(all_S_indicators_d1,2);
inc_probs_S_d2 = mean(all_S_indicators_d2,2);
inc_probs_S_dinf = mean(all_S_indicators_dinf,2);

% Compute subset sizes
subset_size_SS = sum(all_SS_indicators,1);
subset_size_S_d1 = sum(all_S_indicators_d1,1);
subset_size_S_d2 = sum(all_S_indicators_d2,1);
subset_size_S_dinf = sum(all_S_indicators_dinf,1);

% Compute average subset sizes
avg_subset_size_SS = reshape(mean(subset_size_SS,2), [1, length(K_list)]);
avg_subset_size_S_d1 = reshape(mean(subset_size_S_d1,2), [1, length(K_list)]);
avg_subset_size_S_d2 = reshape(mean(subset_size_S_d2,2), [1, length(K_list)]);
avg_subset_size_S_dinf = reshape(mean(subset_size_S_dinf,2), [1, length(K_list)]);

% Error estimates for average subset sizes

z_alpha = 1.96;
se_avg_subset_size_SS = z_alpha*reshape(std(subset_size_SS,0,2), [1, length(K_list)])/sqrt(macroreps);
se_avg_subset_size_S_d1 = z_alpha*reshape(std(subset_size_S_d1,0,2), [1, length(K_list)])/sqrt(macroreps);
se_avg_subset_size_S_d2 = z_alpha*reshape(std(subset_size_S_d2,0,2), [1, length(K_list)])/sqrt(macroreps);
se_avg_subset_size_S_dinf = z_alpha*reshape(std(subset_size_S_dinf,0,2), [1, length(K_list)])/sqrt(macroreps);

%% Compute the true objective function

cost = 3; % per unit order cost 
sell_price = 9; % per unit sale revenue
salvage = 1; % per unit salvage revenue
shortage = 1; % per unit shortage cost
wbl_scale = 50;
wbl_shape = 2; 

feas_region = [1:200]';
true_mean = zeros(length(feas_region),1);
neg_profit = @(D,Q) (cost*Q + shortage*(D - min(D,Q)) - sell_price*min(D,Q) - salvage*(Q - min(D,Q)));
for i = 1:length(true_mean)
    Q = feas_region(i);
    true_mean(i) = integral(@(D) neg_profit(D, Q).*wblpdf(D, wbl_scale, wbl_shape), 0, Inf);
end
[opt_val, opt_index] = min(true_mean);

opt_gap = true_mean - opt_val;
repmat_opt_gap = repmat(opt_gap, [1, macroreps, length(K_list)]);

% Compute optimality gaps
opt_gaps_SS = all_SS_indicators.*repmat_opt_gap;
opt_gaps_S_d1 = all_S_indicators_d1.*repmat_opt_gap;
opt_gaps_S_d2 = all_S_indicators_d2.*repmat_opt_gap;
opt_gaps_S_dinf = all_S_indicators_dinf.*repmat_opt_gap;

% Compute (screening) power -> average optimality gap
% Avoid division by zero when subset size = 0
ss_eps = 0.000001;
power_SS  = reshape(mean(sum(opt_gaps_SS, 1)./(subset_size_SS+eps),2), [1, length(K_list)]);
power_S_d1 = reshape(mean(sum(opt_gaps_S_d1, 1)./(subset_size_S_d1+eps),2), [1, length(K_list)]);
power_S_d2 = reshape(mean(sum(opt_gaps_S_d2, 1)./(subset_size_S_d2+eps),2), [1, length(K_list)]);
power_S_dinf = reshape(mean(sum(opt_gaps_S_dinf, 1)./(subset_size_S_dinf+eps),2), [1, length(K_list)]);

% Compute coverage
coverage_SS = inc_probs_SS(opt_index,:);
coverage_S_d1 = inc_probs_S_d1(opt_index,:);
coverage_S_d2 = inc_probs_S_d2(opt_index,:);
coverage_S_dinf = inc_probs_S_dinf(opt_index,:);

%% Compute S(X) (Fixed K)
exp_set = [20; 60; 100; 140; 180];

addpath('..\..\src');
SX_indicators = construct_det_subset(feas_region, exp_set, true_mean, 'convex', '');
%SX_indicators = construct_det_subset(feas_region, exp_set, true_mean, 'lipschitz_proj', 7);
rmpath('..\..\src');

diff_indicators = diff(SX_indicators);
switch_on = find(diff(SX_indicators) == 1)+1;
switch_off = find(diff(SX_indicators) == -1);

if min(switch_off) < min(switch_on)
    switch_on = [1; switch_on];
end
if max(switch_off) < max(switch_on)
    switch_off = [switch_off; 200];
end

% figure
% hold on
% for i = 1:length(switch_on)
%     x = [switch_on(i)-.5, switch_off(i)+.5, switch_off(i)+.5, switch_on(i)-.5];
%     y = [0, 0, 1, 1]; 
%     p=patch(x,y,'r','LineStyle','none','FaceAlpha',0.5);
% end
% hold off

%% Compute S(X) (Fixed N)

addpath('..\..\src');

SX_indicators = zeros(card_feas_region, length(K_list));

for i = 1:length(K_list)
    K = K_list(i);
    exp_set = round((1:K)'*(200/K) - (200/(2*K)));
    %SX_indicators(:,i) = construct_det_subset(feas_region, exp_set, true_mean, 'convex', '');
    SX_indicators(:,i) = construct_det_subset(feas_region, exp_set, true_mean, 'lipschitz_proj', 7);
end

rmpath('..\..\src');


%% Plotting Setup
string_names = {'STB', 'PO: $d^1$', 'PO: $d^2$', 'PO: $d^{\infty}$'};
colors = {'k-', 'b-', 'g-', 'm-'};
rgb_blue = [0, 0.4470, 0.7410];	
rgb_red = [0.8500, 0.3250, 0.0980]; 
rgb_yellow = [0.9290, 0.6940, 0.1250];
rgb_purple = [0.4940, 0.1840, 0.5560];

%% Plot average sample size

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)
xlim(x_axis_limits)
ylim([0,200])
xlabel(x_axis_label,'interpreter','latex')
ylabel('Average Subset Size ($|\mathcal{S}|$)','interpreter','latex')
%title(plt_title)
hold on
h1 = plot(x_axis_pts, avg_subset_size_SS, 'o-', 'markerfacecolor', rgb_blue, 'MarkerSize', 6, 'LineWidth', 1.5);
h2 = plot(x_axis_pts, avg_subset_size_S_d1, '^-', 'markerfacecolor', rgb_red, 'MarkerSize', 6, 'LineWidth', 1.5);
h3 = plot(x_axis_pts, avg_subset_size_S_d2, 's-', 'markerfacecolor', rgb_yellow, 'MarkerSize', 6, 'LineWidth', 1.5);
h4 = plot(x_axis_pts, avg_subset_size_S_dinf, 'p-', 'markerfacecolor', rgb_purple, 'MarkerSize', 6, 'LineWidth', 1.5);
%line([min(x_axis_pts), max(x_axis_pts)], [sum(SX_indicators), sum(SX_indicators)], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5)
plot(x_axis_pts, sum(SX_indicators,1), 'k:', 'LineWidth', 1.5);
legend([h1, h2, h3, h4], string_names, 'location', 'southeast', 'Interpreter', 'latex');
legend boxoff
box on
hold off

print(['avg_ss_fixedN_ctsnews_',fn_property],'-dpng','-r500')

%x_axis_limits = [min(K_list), max(K_list)];
    %x_axis_pts = K_list;
    %x_axis_label 
%% Plot power

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)
xlim(x_axis_limits)
%ylim([0,200])
xlabel(x_axis_label,'interpreter','latex')
ylabel('Average Average Optimality Gap','interpreter','latex')
title(plt_title)
hold on
plot(x_axis_pts, power_SS, 'o-', 'markerfacecolor', rgb_blue, 'MarkerSize', 6, 'LineWidth', 1.5);
plot(x_axis_pts, power_S_d1, '^-', 'markerfacecolor', rgb_red, 'MarkerSize', 6, 'LineWidth', 1.5);
plot(x_axis_pts, power_S_d2, 's-', 'markerfacecolor', rgb_yellow, 'MarkerSize', 6, 'LineWidth', 1.5);
plot(x_axis_pts, power_S_dinf, 'p-', 'markerfacecolor', rgb_purple, 'MarkerSize', 6, 'LineWidth', 1.5);
legend(string_names, 'location', 'southeast', 'Interpreter', 'latex');
legend boxoff
hold off

%% Plot coverage

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)
xlim(x_axis_limits)
ylim([0,1])
xlabel(x_axis_label,'interpreter','latex')
ylabel('Coverage ($P(x^* \in S)$)','interpreter','latex')
title(plt_title)
hold on
plot(x_axis_pts, coverage_SS, 'o-', 'markerfacecolor', rgb_blue, 'MarkerSize', 6, 'LineWidth', 1.5);
plot(x_axis_pts, coverage_S_d1, '^-', 'markerfacecolor', rgb_red, 'MarkerSize', 6, 'LineWidth', 1.5);
plot(x_axis_pts, coverage_S_d2, 's-', 'markerfacecolor', rgb_yellow, 'MarkerSize', 6, 'LineWidth', 1.5);
plot(x_axis_pts, coverage_S_dinf, 'p-', 'markerfacecolor', rgb_purple, 'MarkerSize', 6, 'LineWidth', 1.5);
legend(string_names, 'location', 'southeast', 'Interpreter', 'latex');
legend boxoff
hold off

%% Plot single 1-D P(x0 in S) for d2

K = 5;
N = 400;
M = 3000;
alpha = 0.05;
load(['ctsnews_N=',num2str(N),'_K=',num2str(K),'_M=',num2str(M),'_iid_',fn_property,'.mat'],'exp_set');

grey_rgb = (192/256)*ones(1,3);
dark_grey_rgb = (128/256)*ones(1,3);
figure
set(gca, 'FontSize', 14, 'LineWidth', 2)

xlim([0,200])
ylim([0,1.005])
xlabel('Solution ($x$)','interpreter','latex')
ylabel('$P(x_0 \in \mathcal{S})$','interpreter','latex')
%title(string_names{1},'interpreter','latex')

hold on
plot([1:200], (true_mean - min(true_mean))./(max(true_mean) - min(true_mean)), 'color', dark_grey_rgb, 'LineWidth', 1)
line([0,200], [1-alpha, 1-alpha], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5)
C = sum(all_S_indicators_d2(:,:,1),2)'/M;
%scatter(feas_region(:,1), C, 'ko','markerfacecolor','k')
plot(feas_region(:,1), C, 'b-', 'LineWidth', 2);
plot(exp_set,zeros(1,size(exp_set,2)),'kx','markerfacecolor','k', 'MarkerSize', 12, 'LineWidth', 1)

for i = 1:length(switch_on)
    x = [switch_on(i)-.5, switch_off(i)+.5, switch_off(i)+.5, switch_on(i)-.5];
    y = [0, 0, 1, 1]; 
    p=patch(x,y,'b','LineStyle','none','FaceAlpha',0.3);
end
box on
hold off

print(['inc_probs_ctsnews_',fn_property,'_K=5_N=400'],'-dpng','-r500')

%% Plot 1-D P(x0 in S)

K = 5;
N = 400;
M = 1;
alpha = 0.05;

load(['ctsnews_N=',num2str(N),'_K=',num2str(K),'_iid_',fn_property,'.mat'],'exp_set');

grey_rgb = (192/256)*ones(1,3);

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)

%
subplot(1,4,1)

axis square
xlim([0,200])
ylim([0,1])
xlabel('Solution ($x$)','interpreter','latex')
ylabel('$P(x_0 \in S)$','interpreter','latex')
title(string_names{1},'interpreter','latex')

hold on
yyaxis left
C = sum(all_SS_indicators(:,:,1),2)'/200;
%scatter(feas_region(:,1), C, 'ko','markerfacecolor','k')
plot(feas_region(:,1), C, 'b-', 'LineWidth', 2)
line([0,200], [1-alpha, 1-alpha], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5)
yyaxis right
plot([1:200],true_mean, 'color', grey_rgb, 'LineWidth',1)
plt = gca;
plt.YAxis(1).Color = 'k';
plt.YAxis(2).Color = 'k';
set(gca,'ytick',[]);
set(gca,'ycolor',[1 1 1])
hold off

%
subplot(1,4,2)

axis square
xlim([0,200])
ylim([0,1])
xlabel('Solution ($x$)','interpreter','latex')
title(string_names{2},'interpreter','latex')

hold on
yyaxis left
C = sum(all_S_indicators_d1(:,:,1),2)'/200;
%scatter(feas_region(:,1), C, 'ko','markerfacecolor','k')
plot(feas_region(:,1), C, 'b-', 'LineWidth', 2)
plot(exp_set,zeros(1,size(exp_set,2)),'kx','markerfacecolor','k', 'MarkerSize', 8, 'LineWidth', 1)
line([0,200], [1-alpha, 1-alpha], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5)
yyaxis right
plot([1:200],true_mean, 'color', grey_rgb, 'LineWidth',1)
plt = gca;
plt.YAxis(1).Color = 'k';
plt.YAxis(2).Color = 'k';
set(gca,'ytick',[]);
set(gca,'ycolor',[1 1 1])
hold off

%
subplot(1,4,3)

axis square
xlim([0,200])
ylim([0,1])
xlabel('Solution ($x$)','interpreter','latex')
title(string_names{3},'interpreter','latex')

hold on
yyaxis left
C = sum(all_S_indicators_d2(:,:,1),2)'/200;
%scatter(feas_region(:,1), C, 'ko','markerfacecolor','k')
plot(feas_region(:,1), C, 'b-', 'LineWidth', 2)
plot(exp_set,zeros(1,size(exp_set,2)),'kx','markerfacecolor','k', 'MarkerSize', 8, 'LineWidth', 1)
line([0,200], [1-alpha, 1-alpha], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5)
yyaxis right
plot([1:200],true_mean, 'color', grey_rgb, 'LineWidth',1)
plt = gca;
plt.YAxis(1).Color = 'k';
plt.YAxis(2).Color = 'k';
set(gca,'ytick',[]);
set(gca,'ycolor',[1 1 1])
hold off

%
subplot(1,4,4)

axis square
xlim([0,200])
ylim([0,1])
xlabel('Solution ($x$)','interpreter','latex')
title(string_names{4},'interpreter','latex')

hold on
yyaxis left
C = sum(all_S_indicators_dinf(:,:,1),2)'/200;
%scatter(feas_region(:,1), C, 'ko','markerfacecolor','k')
plot(feas_region(:,1), C, 'b-', 'LineWidth', 2)
plot(exp_set,zeros(1,size(exp_set,2)),'kx','markerfacecolor','k', 'MarkerSize', 8, 'LineWidth', 1)
line([0,200], [1-alpha, 1-alpha], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5)
yyaxis right
plot([1:200],true_mean, 'color', grey_rgb, 'LineWidth',1)
plt = gca;
plt.YAxis(1).Color = 'k';
plt.YAxis(2).Color = 'k';
set(gca,'ytick',[]);
set(gca,'ycolor',[1 1 1])
hold off

%% Make CRN plots
myVars = {'SS_indicators_CRN', 'S_indicators_dcrn'};
%load(['ctsnews_N=400_K=5_M=3000_crn_convex.mat'],myVars{:});
%load(['ctsnews_N=400_K=5_M=3000_crn_lipschitz_proj.mat'],myVars{:});
load(['ctsnews_N=400_K=5_M=3000_crn_lipschitz_proj_reg.mat'],myVars{:});

% Compute P(x0 in S)
inc_probs_SS_CRN = mean(SS_indicators_CRN,2);
inc_probs_S_dcrn = mean(S_indicators_dcrn,2);

% Compute subset sizes
subset_size_SS_CRN = sum(SS_indicators_CRN,1);
subset_size_S_dcrn = sum(S_indicators_dcrn,1);

% Compute average subset sizes
avg_subset_size_SS_CRN = mean(subset_size_SS_CRN,2);
avg_subset_size_S_dcrn = mean(subset_size_S_dcrn,2);

cost = 3; % per unit order cost 
sell_price = 9; % per unit sale revenue
salvage = 1; % per unit salvage revenue
shortage = 1; % per unit shortage cost
wbl_scale = 50;
wbl_shape = 2; 

feas_region = [1:200]';
true_mean = zeros(length(feas_region),1);
neg_profit = @(D,Q) (cost*Q + shortage*(D - min(D,Q)) - sell_price*min(D,Q) - salvage*(Q - min(D,Q)));
for i = 1:length(true_mean)
    Q = feas_region(i);
    true_mean(i) = integral(@(D) neg_profit(D, Q).*wblpdf(D, wbl_scale, wbl_shape), 0, Inf);
end
[opt_val, opt_index] = min(true_mean);

exp_set = [20; 60; 100; 140; 180];

addpath('..\..\src');
%SX_indicators = construct_det_subset(feas_region, exp_set, true_mean, 'convex', '');
SX_indicators = construct_det_subset(feas_region, exp_set, true_mean, 'lipschitz_proj', 7);
rmpath('..\..\src');

diff_indicators = diff(SX_indicators);
switch_on = find(diff(SX_indicators) == 1)+1;
switch_off = find(diff(SX_indicators) == -1);

if min(switch_off) < min(switch_on)
    switch_on = [1; switch_on];
end
if max(switch_off) < max(switch_on)
    switch_off = [switch_off; 200];
end

string_names = {'STB w/ CRN', 'PO: $d^\mathrm{CRN}$'};

K = 5;
N = 400;
M = 3000;
alpha = 0.05;

%load(['ctsnews_N=',num2str(N),'_K=',num2str(K),'_crn_lipschitz_proj.mat'],'exp_set');

grey_rgb = (192/256)*ones(1,3);
dark_grey_rgb = (128/256)*ones(1,3);

figure
set(gca, 'FontSize', 14, 'LineWidth', 2)

%

axis square
xlim([0,200])
ylim([0,1.005])
xlabel('Solution ($x$)','interpreter','latex')
ylabel('$P(x_0 \in \mathcal{S})$','interpreter','latex')
%title(string_names{1},'interpreter','latex')
% 
% hold on
% yyaxis left
% C = sum(SS_indicators_CRN,2)'/200;
% %C = sum(S_indicators_dcrn,2)'/200;
% %scatter(feas_region(:,1), C, 'ko','markerfacecolor','k')
% plot(feas_region(:,1), C, 'b-', 'LineWidth', 2)
% line([0,200], [1-alpha, 1-alpha], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5)
% yyaxis right
% plot([1:200],true_mean, 'color', grey_rgb, 'LineWidth',1)
% plt = gca;
% plt.YAxis(1).Color = 'k';
% plt.YAxis(2).Color = 'k';
% set(gca,'ytick',[]);
% set(gca,'ycolor',[1 1 1])
% hold off

hold on
plot([1:200], (true_mean - min(true_mean))./(max(true_mean) - min(true_mean)), 'color', dark_grey_rgb, 'LineWidth', 1)
line([0,200], [1-alpha, 1-alpha], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5)
%C = sum(SS_indicators_CRN,2)'/M;
C = sum(S_indicators_dcrn,2)'/M;
plot(feas_region(:,1), C, 'b-', 'LineWidth', 2);
plot(exp_set,zeros(1,size(exp_set,2)),'kx','markerfacecolor','k', 'MarkerSize', 12, 'LineWidth', 1)

for i = 1:length(switch_on)
    x = [switch_on(i)-.5, switch_off(i)+.5, switch_off(i)+.5, switch_on(i)-.5];
    y = [0, 0, 1, 1]; 
    p=patch(x,y,'b','LineStyle','none','FaceAlpha',0.3);
end
box on
hold off

%print(['inc_probs_ctsnews_lipschitz_crn_STB_K=5_N=4000'],'-dpng','-r500')


% %
% figure
% 
% axis square
% xlim([0,200])
% ylim([0,1])
% xlabel('Solution ($x$)','interpreter','latex')
% title(string_names{2},'interpreter','latex')
% 
% hold on
% yyaxis left
% C = sum(S_indicators_dcrn,2)'/200;
% %scatter(feas_region(:,1), C, 'ko','markerfacecolor','k')
% plot(feas_region(:,1), C, 'b-', 'LineWidth', 2)
% plot(exp_set,zeros(1,size(exp_set,2)),'kx','markerfacecolor','k', 'MarkerSize', 8, 'LineWidth', 1)
% line([0,200], [1-alpha, 1-alpha], 'Color', 'black', 'LineStyle', ':', 'LineWidth', 1.5)
% yyaxis right
% plot([1:200],true_mean, 'color', grey_rgb, 'LineWidth',1)
% plt = gca;
% plt.YAxis(1).Color = 'k';
% plt.YAxis(2).Color = 'k';
% set(gca,'ytick',[]);
% set(gca,'ycolor',[1 1 1])
% hold off
