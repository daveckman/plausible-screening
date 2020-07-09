function [outputs] = sS_2prod_oracle(oracle_rngs, solution, n_reps)
% Take n_reps replications at a given solution.
% Return the outputs as a row vector.
% Input solution is a row vector.

% sS_oracle evaluates a 30-day period of an (s,S) inventory reorder policy
% for two products. Each products parameters are similar to those in Koenig
% and Law (1985). See Pei and Nelson research for more details.

% Parameters
t = 30; % t is the time horizon
o = 32; % o is the order cost
u = 3; % u is the unit shipping cost
h = 1; % h is the holding cost, per unit, per period
b = 5; % b is the back-order cost, per unit, per period
lambda = 25; %lambda is the mean of Poisson-distributed demand

% Initialize for storage
outputs = zeros(1,n_reps);
avg_cost = zeros(2, n_reps);

for system = 1:2
    % Unpack (s,S) policy to evaluate
    ss = solution(1 + 2*(system-1)); % s
    QQ = solution(2 + 2*(system-1)); % S - s
    
    % Unpack random number stream
    demand_stream = oracle_rngs{system};
    RandStream.setGlobalStream(demand_stream);
    
    for j = 1:n_reps

        % Initialize
        I = ss + QQ;
        C = 0;
        
        % Generate all monthly demands up front
        Ds = poissrnd(lambda,[t,1]);
        
        % Simulate a t-month period
        
        for k = 1:t
            
            reorder = (I < ss);
            quantity = ss + QQ - I;
            I = I + reorder * quantity;
            C = C + reorder * (o + u * quantity);
            I = I - Ds(k);
            multiplier = h * (I > 0) - b * (I < 0);
            C = C + I * multiplier;
            
        end
        
        % Calculate average monthly cost
        avg_cost(system, j) = C/t;
        
    end
end

for j = 1:n_reps
    outputs(j) = avg_cost(1, j) + avg_cost(2, j) - sqrt((solution(3) - 18)^2 + (solution(4) - 35)^2); 
end

end

