function [outputs] = newsvendor_oracle(oracle_rngs, solution, n_reps)
% Take n_reps replications at a given solution.
% Return the outputs as a row vector.
% Input solution is a row vector.

% newsvendor_oracle simulates a day's demand for the classical newsvendor
% problem (one-dimensional). It returns the loss (negative profit) for a given order
% quantity. Demand follows a Poisson distribution.

% Unpack random number streams
demand_stream = oracle_rngs{1};
RandStream.setGlobalStream(demand_stream);

% Unpack order quantity to evaluate
order_quantity = solution;

% Initialize for storage
outputs  = zeros(1,n_reps);

% Parameters
cost = 3; % per unit order cost 
sell_price = 9; % per unit sale revenue
salvage = 1; % per unit salvage revenue
lambda = 50; % average daily demand

% Generate outputs
for j = 1:n_reps
    
    % New substream for each replication
    demand_stream.Substream = j;
    
    % Generate daily demand
    demand = poissrnd(lambda);
    
    % Calculate loss (negative profit)
    order_cost = cost * order_quantity;
    sales_revenue = sell_price * min(demand, order_quantity);
    salvage_revenue = salvage * (order_quantity - min(demand, order_quantity));
    outputs(j) = order_cost - sales_revenue - salvage_revenue;
    
end