function plot_sample_size_ecdfs(all_S_indicators, string_names, colors)
% Plot figure with sample sizes ecdf curves

num_algs = size(all_S_indicators, 2);
ecdfs = cell(1, num_algs);

card_feas_region = size(all_S_indicators{1}, 1);

for i = 1:num_algs
    subsize = sum(all_S_indicators{i}, 1);
    ecdf = zeros(1, card_feas_region + 1);
    for j = 0:card_feas_region
        ecdf(j+1) = mean(subsize <= j);
    end
    ecdfs{i} = ecdf;
end

%figure
figure('Position', [0, 0, 800, 400])

set(gca, 'FontSize', 14, 'LineWidth', 2)

xlim([0,card_feas_region])
ylim([0,1])
xlabel('Subset Size ($s$)', 'Interpreter', 'latex')
ylabel('P($|S| \leq s$)', 'Interpreter', 'latex')

hold on
for i = 1:num_algs
    stairs([0:card_feas_region], ecdfs{i}, colors{i}, 'LineWidth', 2);
end
hold off

leg1 = legend(string_names, 'location', 'northwest', 'Interpreter', 'latex');
legend boxoff
set(leg1,'Interpreter','latex');
set(gca, 'FontSize', 14, 'LineWidth', 2)
%set(gcf,'units','inches','position',[1,1,9,4.5])

hold off