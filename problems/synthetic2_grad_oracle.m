function [outputs, gradients] = synthetic2_grad_oracle(oracle_rngs, solution, n_reps)
% Take n_reps replications at a given solution.
% Return the outputs as a row vector. (1 x n)
% Return the gradients as column vectors. (d x n)
% Input solution is a row vector.

% Simulate noisy observations and gradients from
% mu(x) = x^2 - x - x*y - y + y^2 + 1
% GRAD mu(x) = [2x - y - 1; 2y - x - 1] 

% Unpack random number streams
noise_stream = oracle_rngs{1};
RandStream.setGlobalStream(noise_stream);

% Unpack solution
x = solution(1);
y = solution(2);

% Compute magnitude of local Gaussian noise
mu_noise = 0.5 + sqrt((x-1)^2 + (y-1)^2);
grad_noise = 0.5 + sqrt((x-1)^2 + (y-1)^2);
rho_mu_x = 0.5;
rho_mu_y = 0.3;
rho_x_y = -0.2;

% (Y, G) will be generated according to MVN(MU, SIGMA)
MU = [x^2 - x - x*y - y + y^2 + 1, 2*x - y - 1, 2*y - x - 1];
SIGMA = [mu_noise, rho_mu_x*sqrt(mu_noise)*sqrt(grad_noise), rho_mu_y*sqrt(mu_noise)*sqrt(grad_noise);
    rho_mu_x*sqrt(mu_noise)*sqrt(grad_noise), grad_noise, rho_x_y*sqrt(grad_noise)*sqrt(grad_noise);
    rho_mu_y*sqrt(mu_noise)*sqrt(grad_noise), rho_x_y*sqrt(grad_noise)*sqrt(grad_noise), grad_noise];

% Initialize for storage
outputs  = zeros(1,n_reps);
d = length(solution);
gradients = zeros(d,n_reps);

% Generate outputs
for j = 1:n_reps
    
    % New substream for each replication
    noise_stream.Substream = j;
    
    % Generate Z = (Y, G) ~ MVN((mu(x), GRAD mu(x)), 
    Z = mvnrnd(MU, SIGMA);
    outputs(j) = Z(1);
    gradients(:,j) = Z(2:(d+1))';
        
end