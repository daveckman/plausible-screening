clear
clc;

K_list = [5, 10, 20, 50];

myVars = {'SS_indicators', 'S_indicators_d1', 'S_poly_indicators_d1', 'S_indicators_d2', 'S_poly_indicators_d2', 'S_indicators_dinf', 'S_poly_indicators_dinf'};
string_names = {'ExtSTB', 'PO: $d^1$', 'PO (relaxed): $d^1$', 'PO: $d^2$', 'PO (relaxed): $d^2$', 'PO: $d^{\infty}$', 'PO (relaxed): $d^{\infty}$'};
colors = {'k:', 'b-', 'b:', 'g-', 'g:', 'm-', 'm:'};
    
for i = 1:length(K_list)
    load(['ctsnews_N=400_K=',num2str(K_list(i)),'_iid_lipschitz.mat'],myVars{:});
    all_S_indicators = {SS_indicators, S_indicators_d1, S_poly_indicators_d1, S_indicators_d2, S_poly_indicators_d2, S_indicators_dinf, S_poly_indicators_dinf};
    
    plot_sample_size_ecdfs(all_S_indicators, string_names, colors)
end

