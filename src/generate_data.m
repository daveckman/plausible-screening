function [sample_mean, sample_var, sample_pair_var] = generate_data(m, oracle_string, oracle_n_rngs, exp_set, n_vec, discrep_string)
% Generate data from solutions in the experiment set.
% Take n_i replications at solution x_i for i = 1, ... k where k is the
% cardinality of the experimental set X.
% m = macroreplication number

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
        oracle_rngs{r} = RandStream.create('mrg32k3a', 'NumStreams', m*oracle_n_rngs, 'StreamIndices', (m - 1)*oracle_n_rngs + r);
    end
    
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
        
    end % end for
    
    % Calculate summary statistics
    sample_mean = mean(outputs,2);
    sample_var = cov(outputs');
    
    % !! TO AVOID SINGULARITIES !!
    sample_var = sample_var + 1e-6*eye(k);
    
    % Calculate pairwise difference variances (for ESTB with CRN)
    sample_pair_var = zeros(k, k);
    for j = 1:k
        for l = 1:k
            sample_pair_var(j, l) = 1/(n_vec(1) - 1)*sum((outputs(j,:) - outputs(l,:) - (sample_mean(j) - sample_mean(l))).^2);
        end
    end
    
else % otherwise do independent sampling
    
    % Initialize for storage
    sample_mean = zeros(k,1);
    sample_var = zeros(k,1);
    
    for i = 1:k
        
        % Extract solution x_i and sample size n_i
        x_i = exp_set(i,:);
        n_i = n_vec(i);

        % Set up distinct random number streams to use
        oracle_rngs = cell(1, oracle_n_rngs);
        for r = 1:oracle_n_rngs
            oracle_rngs{r} = RandStream.create('mrg32k3a', 'NumStreams', m*oracle_n_rngs*k, 'StreamIndices', (m - 1)*oracle_n_rngs*k + oracle_n_rngs*(i - 1) + r);
        end
        
        % Take n_i replications at x_i
        outputs = oracle_handle(oracle_rngs, x_i, n_i);

        % Calculate summary statistics
        sample_mean(i) = mean(outputs);
        sample_var(i) = var(outputs);
                
    end % end for
    
    sample_pair_var = [];
    
end % end if

% Avoid sample variances of zero;
%sample_var(sample_var == 0) = 0.00001;

end

