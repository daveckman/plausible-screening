function [sample_mean, sample_var] = generate_data(oracle_string, oracle_n_rngs, exp_set, n_vec, discrep_string)
% Generate data from solutions in the experiment set.
% Take n_i replications at solution x_i for i = 1, ... k where k is the
% cardinality of the experimental set X.

oracle_handle = str2func(oracle_string);

% Each row of exp_set corresponds to a solution 
k = size(exp_set, 1); % Number of evaluated solutions

if strcmp(discrep_string, 'CRN') == 1
    % Use CRN if CRN discrepancy is specified
    
    % Initialize for storage
    outputs = zeros(k, n_vec(1));
    
    % Create one set of streams to reuse
    oracle_rngs = cell(1, oracle_n_rngs);
    for r = 1:oracle_n_rngs
        oracle_rngs{r} = RandStream.create('mrg32k3a', 'NumStreams', oracle_n_rngs, 'StreamIndices', r);
    end
    
    parfor_progress(k);
    for i = 1:k  
        % Extract solution x_i and sample size n_i
        x_i = exp_set(i,:);
        n_i = n_vec(i);
        
        % Reset each stream to substream 1
        for r = 1:oracle_n_rngs
            oracle_rng = oracle_rngs{r};
            oracle_rng.Substream = 1;
        end
        
        % Take n_i replications at x_i
        outputs(i,:) = oracle_handle(oracle_rngs, x_i, n_i);
        
        parfor_progress;
    end
    parfor_progress(0);
    
    % Calculate summary statistics
    sample_mean = mean(outputs,2);
    sample_var = cov(outputs');
    
else % otherwise do independent sampling
    
    % Initialize for storage
    sample_mean = zeros(k,1);
    sample_var = zeros(k,1);
    
    parfor_progress(k);
    for i = 1:k
        
        % Extract solution x_i and sample size n_i
        x_i = exp_set(i,:);
        n_i = n_vec(i);

        % Set up distinct random number streams to use
        oracle_rngs = cell(1, oracle_n_rngs);
        for r = 1:oracle_n_rngs
            oracle_rngs{r} = RandStream.create('mrg32k3a', 'NumStreams', oracle_n_rngs*k, 'StreamIndices', oracle_n_rngs*(i - 1) + r);
        end
        
        % Take n_i replications at x_i
        outputs = oracle_handle(oracle_rngs, x_i, n_i);
        
        % Calculate summary statistics
        sample_mean(i) = mean(outputs);
        sample_var(i) = var(outputs);
        
        parfor_progress;
        
    end
    parfor_progress(0);
    
end % end if

end

