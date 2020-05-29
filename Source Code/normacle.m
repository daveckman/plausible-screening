function [outputs] = normacle(oracle_rngs, solution, n_reps)
% Take n_reps replications at a given solution.
% Return the outputs as a row vector.
% Input solution is a row vector.

% Normacle evaluates the Euclidean norm of solution with Gaussian error
% Variance of the noise is proportional to the norm of the solution.

% Unpack random number streams
noise_stream = oracle_rngs{1};


% Initialize for storage
outputs = zeros(1, n_reps);

% Generate outputs
for j = 1:n_reps
    
    % New substream for each iteration
    noise_stream.Substream = j;
    
    % Generate noise
    RandStream.setGlobalStream(noise_stream)
    noise = normrnd(0, norm(solution, 2));
    
    % Add noise to expected value function
    outputs(j) = norm(solution, 2) + noise;
end

end

