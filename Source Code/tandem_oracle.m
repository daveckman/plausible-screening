function [outputs] = tandem_oracle(oracle_rngs, solution, n_reps)
% Take n_reps replications at a given solution.
% Return the outputs as a row vector.
% Input solution is a row vector.

% tandem_oracle simulates a tandem production line with manufacturing
% blocking (buffers). It returns the completion time of the 100th product
% given parameters for the service time distributions of each machine.
% See Plambeck et al. (1996) and Shanthikumar and Yao (1989) for details.

% Unpack random number streams
service_stream = oracle_rngs{1};
RandStream.setGlobalStream(service_stream);
% Because a fixed number of products move through the system and are
% serviced according to a FIFO discipline, each product's vector of service
% times at all machines can be generated in advance.

% Unpack configuration of machines' mean cycle times to evaluate
mean_cycle_times = solution;

% Initialize for storage
outputs  = zeros(1,n_reps);

% Description of the tandem production line
%   - tandem configuration of single-server machines
%   - machine cycle (service) times are exponentially distributed
%   - between pairs of machines there are buffers of fixed size
%   - buffer at a machine include the spot for the machine being processed
%   - if buffer is full, upstream machine can become blocked
%   - if buffer is empty, downstream machine can become starved
%   - no arrival process -> first machine has all products available to process

% Other setup (needs to match setup in test_tandem_production.m)
buffers = [Inf, 3]; % vector of buffers in front of each machine
num_machines = size(buffers, 2); % number of machines
num_prod = 100; % number of products to process

% Generate outputs
for j = 1:n_reps
    
    % Let S_n^i be the service time of product n at machine i
    % Let D_n^i be the completion time of product n at machine i

    % Recursion from Shanthikumar and Yao (1989; Eqn (7)):
    %   D_n^i = max([max(D_n^{i-1}, D_{n-1}^i) + S_n^i], D_{n-{b_i}+1}^{i+1})
    %   with D_k^i := 0 for all k <= 0 and D_k^{m+1} := 0 for all k.

    % New substream for each replication
    service_stream.Substream = j;
    
    % Generate all service times up front and store in matrix
    S = zeros(num_prod, num_machines);
    for col = 1:num_machines
        S(:,col) = exprnd(mean_cycle_times(col), [num_prod, 1]);
    end
    
    % Compute completion times (via above recursion) and store in matrix
    D = zeros(num_prod, num_machines);
    D(1,:) = cumsum(S(1,:)); % first product is processed without delays
    for prod = 2:num_prod
        for machine = 1:num_machines
            
            if machine == 1
                comp_prev_machine = 0;
            else
                comp_prev_machine = D(prod, machine - 1);
            end
                                    
            if machine == num_machines
                comp_next_buffer = 0;
            elseif prod - buffers(machine + 1) <= 0
                comp_next_buffer = 0;
            else
                comp_next_buffer = D(prod - buffers(machine + 1), machine + 1); 
            end
            
            D(prod, machine) = max(max(comp_prev_machine, D(prod - 1, machine)) + S(prod, machine), comp_next_buffer);
        end
    end 
    
    % Report completion time of last product
    outputs(j) = D(num_prod, num_machines);
    
end