function plot_1d_subset_prob(feas_region, exp_set, S_indicators, true_mean, xaxlabel, yyaxlabel, title_str, alpha)
% Plot subplot of inclusion probabilities for a given screening algorithm
% Superimpose true objective function

set(gca, 'FontSize', 14, 'LineWidth', 2)

axis square
xlim([min(feas_region), max(feas_region)])
ylim([0,1])
xlabel(xaxlabel, 'interpreter', 'latex')
ylabel('$P(x_0 \in S)$','interpreter','latex')
title(title_str, 'interpreter', 'latex')

hold on
yyaxis left
C = sum(S_indicators,2)'/size(S_indicators,2);
%scatter(feas_region(:,1), C, 'ko','markerfacecolor','k')
plot(feas_region(:,1), C, 'g:', 'LineWidth', 2)
plot(exp_set(:,1),zeros(1,size(exp_set,2)),'kx','markerfacecolor','k', 'MarkerSize', 12, 'LineWidth', 2)
line([min(feas_region), max(feas_region)], [1-alpha, 1-alpha], 'Color', 'black', 'LineStyle', '-', 'LineWidth', 2)
yyaxis right
plot(feas_region(:,1),true_mean,'b--', 'LineWidth',2)
plt = gca;
plt.YAxis(1).Color = 'k';
plt.YAxis(2).Color = 'k';
ylabel(yyaxlabel,'interpreter','latex')
hold off

end

