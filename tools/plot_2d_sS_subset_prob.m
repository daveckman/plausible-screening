function plot_2d_sS_subset_prob(feas_region, exp_set, S_indicators, title_str)
% Make a 2d plot of the inclusion probabilities for the (s, S) problem

set(gca, 'FontSize', 14, 'LineWidth', 2)
xlim([min(feas_region(:,1))-5,max(feas_region(:,1))+5])
ylim([min(feas_region(:,2))-5,max(feas_region(:,2))+5])
xlabel('Reorder quantity ($s$)','interpreter','latex')
ylabel('Order-to quantity ($S$)','interpreter','latex')
title(title_str, 'interpreter', 'latex')

hold on
C = sum(S_indicators,2)'/size(S_indicators,2);
for k = 1:size(feas_region,1)
    rectangle('Position',[feas_region(k,1)-C(k)/2  feas_region(k,2)-C(k)/2 C(k) C(k)],'facecolor',[0.2 0.8 0.8],'linestyle','none','Curvature',1)
end

plot(exp_set(:,1),exp_set(:,2),'kx','markerfacecolor','k', 'MarkerSize', 6, 'LineWidth', 2)
vch = convhull(feas_region(:,1),feas_region(:,2));
plot(feas_region(vch,1),feas_region(vch,2),'k-')

hold off

end

