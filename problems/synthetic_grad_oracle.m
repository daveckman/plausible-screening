function [outputs, gradients] = synthetic_grad_oracle(oracle_rngs, solution, n_reps)
% Take n_reps replications at a given solution.
% Return the outputs as a row vector. (1 x n)
% Return the gradients as column vectors. (d x n)
% Input solution is a row vector.

% Simulate noisy observations and gradients from
% mu(x) = 4*(x-0.5)^2
% GRAD mu(x) = 8*x - 4

% Unpack random number streams
noise_stream = oracle_rngs{1};
RandStream.setGlobalStream(noise_stream);

% Unpack solution
x = solution;

% Compute magnitude of local Gaussian noise
mu_noise = 0.5 + (x-0.5)^2;
grad_noise = 0.5 + (x-0.5)^2;
rho = 0.5;

% (Y, G) will be generated according to MVN(MU, SIGMA)
MU = [4*(x-0.5)^2, 8*x-4];
SIGMA = [mu_noise, rho*sqrt(mu_noise)*sqrt(grad_noise); rho*sqrt(mu_noise)*sqrt(grad_noise), grad_noise];

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
    gradients(:,j) = Z(2);
        
end