function SGX_indicators = construct_det_grad_only_subset(feas_region, exp_set, true_mean, true_grad, opttol) %, fn_props, prop_params)
% Construct the deterministic subset SG(X) for known gradient values.
% No known objective values. Using relaxed Z(x0).

% true_grad is a (|feasregion| x d) matrix, i.e., gradients are row vectors

% FOR CONVEX CASE

%PO_info_handle = str2func(strcat('setup_',fn_props));
k = size(exp_set, 1);
card_feas_region = size(feas_region, 1);
SGX_indicators = zeros(card_feas_region, 1);

[~, exp_set_indices] = ismembertol(exp_set, feas_region, 'ByRows', true); 
true_mean_exp_set = true_mean(exp_set_indices,:);
true_grad_exp_set = true_grad(exp_set_indices,:); 


% SG0(X) = {x0 in feasregion:
%          s.t. max_i {-(x_i - x_0)^T*grad(x_i)} <= 0}

for l = 1:card_feas_region
    x0 = feas_region(l,:);
    
    leftmax = max(-sum((exp_set - repmat(x0, k, 1)).*true_grad_exp_set,2));

    SGX_indicators(l) = (leftmax <= opttol); 
end