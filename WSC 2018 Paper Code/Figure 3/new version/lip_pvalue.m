% Create p-value heat map for Lipschitz constant

% Load in raw data from M=100 macroreplications
load('Data_12_02_2019')

% Parameters
%alpha = 0.05;
%beta = (1 - (1 - alpha)^(1./K))/2;
%t_crit = tinv(1-beta, num_sim - 1);

c_vec = 0:0.05:5;

% Vector of p-values
p_c = zeros(M, length(c_vec));

for m = 1:M

    % Construct LHS constraint matrix for Lipschitz
    A_LHS = zeros(2*K*(K-1), K);
    for i = 1:K
        for j = 1:K
            if i ~= j            
                if j < i
                    row = (i-1)*(K-1) + j;
                elseif j > i
                    row = (i-1)*(K-1) + j - 1;
                end
                A_LHS(row, i) = 1;
                A_LHS(row + K*(K-1), i) = -1;
                A_LHS(row, j) = -1; 
                A_LHS(row + K*(K-1), j) = 1;
            end
        end
    end

    % Just look at the first experiment (m = 1);
    muhat = muhat_by_macro(m,:);
    sigma2hat = sigma2hat_by_macro(m,:);

    % Preliminary matrix multiplications
    A_muhat = A_LHS*muhat';
    A_sigma = abs(A_LHS)*sqrt(sigma2hat'/num_sim);

    for c_index = 1:length(c_vec)

        %disp(c_index)

        c = c_vec(c_index);

        % Construct RHS vector b
        b_RHS = zeros(2*K*(K-1),1);

        for i = 1:K
            for j = 1:K
                if i ~= j            
                    if j < i
                        row = (i-1)*(K-1) + j;
                    elseif j > i
                        row = (i-1)*(K-1) + j - 1;
                    end
                    b_RHS(row) = c*norm(X(i,:) - X(j,:));
                    b_RHS(row + K*(K-1)) = c*norm(X(i,:) - X(j,:));
                end
            end
        end

        % Find min alpha value
        min_beta = tcdf((b_RHS - A_muhat)./A_sigma, num_sim - 1);
        beta_bar = min(min_beta);
        p_c(m, c_index) = min(1, 1 - (1 - 2*beta_bar)^K); % since beta can only be in (0, 1/2);

    end
end
%%

plot(c_vec, p_c);
xlabel('Lipschitz Constant')
ylabel('pseudo p-value')
title('Hypothesis testing for Lipschitz')

%print_correctly('sS_numerical')

