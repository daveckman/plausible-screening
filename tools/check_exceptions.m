function check_exceptions(discrep_string, fn_props, n_vec)
% Check for exceptions with problem specifications

% Check discrepancy string
accept_discreps = {'ell1', 'ell2', 'ellinf', 'CRN'};
if ~any(strcmp(discrep_string, accept_discreps))
    fprintf('\nERROR: "%s" is not a valid standardized discrepancy.\n', discrep_string)
    fprintf('Please specify a valid standardized discrepancy:')
    fprintf('\t%s', accept_discreps{:})
    fprintf('.\n')
    return
end

% Check functional property string
accept_fn_props = {'lipschitz', 'lipschitz_proj', 'convex', 'convex_nearopt'};
if ~any(strcmp(fn_props, accept_fn_props))
    fprintf('\nERROR: "%s" is not a valid functional property.\n', fn_props)
    fprintf('Please specify a valid functional property:')
    fprintf('\t%s', accept_fn_props{:})
    fprintf('.\n')
    return
end

% Check sample size vector for common random numbers
if strcmp(discrep_string, 'CRN') == 1
    % Check that all values in n_vec vector are equal
    if min(n_vec) ~= max(n_vec)
        fprintf('\nERROR: All sample sizes must be equal when using CRN.\n')
        return
    end
    
    % Check if n > k so that sample_var is non-singular
    if n_vec(1) <= length(n_vec)
        fprintf('\nERROR: Common sample size n = %d must be strictly greater than k = %d.\n', n_vec(1), length(n_vec))
        return
    end
end


end

