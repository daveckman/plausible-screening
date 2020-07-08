function [outputs] = sS_oracle(oracle_rngs, solution, n_reps)
% Take n_reps replications at a given solution.
% Return the outputs as a row vector.
% Input solution is a row vector.

% sS_oracle evaluates a 100-day period of an (s,S) inventory reorder policy
% See Koenig and Law (1985) for the original problem and Plumlee and Nelson
% (2018) for the details.

% Unpack random number streams
demand_stream = oracle_rngs{1};
RandStream.setGlobalStream(demand_stream);

% Unpack (s,S) policy to evaluate
s = solution(1);
S = solution(2);

% Initialize for storage
outputs  = zeros(1,n_reps);

% Generate outputs
for j = 1:n_reps
    
    % New substream for each replication
    demand_stream.Substream = j;
    
    % Initialize
    I = S;
    cost = 0;
    
    % Generate all monthly demands up front
    Ds = poissrnd(25,[100, 1]);
    
    % Simulate a 100-month period
    for k = 1:100
        J = I;
        D = Ds(k);
        if J <= s % If stockout... reorder.
            cost = cost + 32 + 3*(S-J);
            J = S;
        end
        if J >= D % If inventory... incur per-unit holding cost.
            cost = cost + (J-D);
        else % If no inventory... incur per-unit shortage cost.
            cost = cost + 5*(D-J);
        end
        I = J-D;
    end
    
    % Calculate average monthly cost
    outputs(j) = cost/100;
    
end