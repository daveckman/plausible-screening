function [outputs, gradients] = cts_newsvendor_grad_oracle(oracle_rngs, solution, n_reps)
% Take n_reps replications at a given solution.
% Return the outputs as a row vector. (1 x n)
% Return the gradients as column vectors. (d x n)
% Input solution is a row vector.

% cst_newsvendor_oracle simulates a day's demand for the classical newsvendor
% problem (one-dimensional). It returns the loss (negative profit) for a given order
% quantity and its gradient. Demand follows a Poisson distribution.

% Unpack random number streams
demand_stream = oracle_rngs{1};
RandStream.setGlobalStream(demand_stream);

% Unpack order quantity to evaluate
order_quantity = solution;

% Initialize for storage
outputs  = zeros(1,n_reps);
d = length(solution);
gradients = zeros(d,n_reps);

% Parameters
cost = 3; % per unit order cost 
sell_price = 9; % per unit sale revenue
salvage = 1; % per unit salvage revenue
shortage = 1; % per unit shortage cost
wbl_scale = 50; % shape parameter of Weibull distribution for daily demand
wbl_shape = 2; % shape parameter of Weibull distribution for daily demand
% avg daily demand = wbl_scale*gamma(1 + 1/wbl_shape)

% Generate outputs
for j = 1:n_reps
    
    % New substream for each replication
    demand_stream.Substream = j;
    
    % Generate daily demand
    demand = wblrnd(wbl_scale, wbl_shape);
    
    % Calculate loss (negative profit)
    order_cost = cost * order_quantity;
    sales_revenue = sell_price * min(demand, order_quantity);
    salvage_revenue = salvage * (order_quantity - min(demand, order_quantity));
    shortage_cost = shortage * (demand - min(demand, order_quantity));
    outputs(j) = order_cost + shortage_cost - sales_revenue - salvage_revenue;
    
    gradients(:,j) = cost - sell_price*(demand > order_quantity) - salvage*(demand < order_quantity) - shortage*(demand > order_quantity);
    
end