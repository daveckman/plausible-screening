% tandem production line problem
oracle_string = 'tandem_budget_oracle';
oracle_n_rngs = 1;

% Allocate budget = 50 resources across 5 machines in integral amounts
budget = 50;
num_machines = 5;

% # feasible solutions = (5 multichoose 50) = 316,251
% Because of tight budget constraint (equality), reduce the dim to d = 4.
feas_region = zeros(nchoosek(num_machines + budget -1, budget), num_machines - 1);
row = 1;
for i1 = 0:budget
    for i2 = 0:(budget - i1)
        for i3 = 0:(budget - i1 - i2)
            for i4 = 0:(budget - i1 - i2 - i3)
                feas_region(row,:) = [i1, i2, i3, i4];
                row = row + 1;
            end
        end
    end
end

% For evaluated solutions, allocate budget = 50 across machines in
% multiples of 10 units
frac_budget = 5;
multiple = 10;
% # evaluated solutions = (5 multichoose 5) = 126 = .04% of feas solns

exp_set = zeros(nchoosek(num_machines + frac_budget - 1, frac_budget),4);
row = 1;
for i1 = 0:frac_budget
    for i2 = 0:(frac_budget - i1)
        for i3 = 0:(frac_budget - i1 - i2)
            for i4 = 0:(frac_budget - i1 - i2 - i3)
                exp_set(row,:) = multiple*[i1, i2, i3, i4];
                row = row + 1;
            end
        end
    end
end

% !!! for faster run time
feas_region = exp_set;

k = size(exp_set, 1);

n_vec = 10*ones(k, 1); % col vector
alpha = 0.05; % Confidence level = 1-alpha
discrep_string = 'ell1'; % {'ell1', 'ell2', 'ellinf', 'CRN'}
fn_props = 'convex'; % {'convex', 'lipschitz', 'lipschitz_proj}
prop_params = []; % gamma for Lipschitz constant
LP_solver_string = 'glpk'; % {'MATLAB', 'glpk'}