function D_x0 = calc_min_std_discrep_gradinf(A, C, b, sample_mean, sample_mean_grad, sample_full_cov, n_vec)

% 1. Take polyhedral representation of P.
% 2. Formulate a mathematical program for minimizing the standardized discrepancy.
% 3. Solve the mathematical program and return the minimum standardized discrepancy.

% Determine dimensions of polyhedral representation of M(x_0) = proj_m(P)
% where P = {(m, w): A*m + C*w <= b}
p = size(A,1); % Number of constraints
q = size(C,2); % Number of unprojected components
k = size(sample_mean_grad,1); % Number of solutions in experimental set
d = size(sample_mean_grad,2); % Dimension of solution space

% Concatenate sample means and gradients
concat = [sample_mean'; sample_mean_grad'];
zetahat = reshape(concat, k*(d+1), 1);

%discrep(z, zetahat, Psihat) = max_{i=1,...,k} n_i*(zetahat_i -
%z_i)^T*inv(Psihat_i)*(zetahata_i - z_i)
g = @(z, i) n_vec(i)*(zetahat(((i-1)*(d+1)+1):(i*(d+1))) - z(((i-1)*(d+1)+1):(i*(d+1))))'*inv(sample_full_cov(:,:,i))*(zetahat(((i-1)*(d+1)+1):(i*(d+1))) - z(((i-1)*(d+1)+1):(i*(d+1))));
discrep = @(z) max([g(z,1), g(z,2), g(z,3), g(z,4), g(z,5)]);

% initial solution
z0 = zeros(k*(d+1)+1,1);

% constraints
A_convp = [A, C];
b_convp = b;

% Compute minimum standardized discrepancy
% Minimize max of quadratic functions
options = optimoptions('fmincon','Display','none','MaxIter',500);
[~, D_x0] = fmincon(discrep, z0, A_convp, b_convp, [], [], [], [], [], options);

end