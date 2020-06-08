% For a single macroreplication, plot a curve showing the number of
% solutions in S if a Lipschitz constant of c is used (as fn of c).

tic;

% Load in raw data from M=100 macroreplications
load('Data_12_02_2019')

c = 1; % will be cancelled out later
alpha = 0.05;
beta = (1 - (1 - alpha)^(1./K))/2;
t_crit = tinv(1-beta, num_sim - 1);

% Vector for storing tight c values
c_screen = zeros(1, length(scrX));

% Just look at the first experiment
muhat = muhat_by_macro(1,:);
sigma2hat = sigma2hat_by_macro(1,:);

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

% Preliminary matrix multiplications
A_muhat = A_LHS*muhat';
A_sigma = abs(A_LHS)*sqrt(sigma2hat'/num_sim);

for x0_index = 1:length(scrX)

    % Construct RHS vector
    b_RHS = zeros(2*K*(K-1),1);
    for i = 1:K
        for j = 1:K
            if i ~= j            
                if j < i
                    row = (i-1)*(K-1) + j;
                elseif j > i
                    row = (i-1)*(K-1) + j - 1;
                end
                % c term that was multiplied out front will be cancelled
                b_RHS(row) = c*min(norm(X(i,:) - X(j,:)), norm(X(i,:) - scrX(x0_index,:)));
                b_RHS(row + K*(K-1)) = c*min(norm(X(i,:) - X(j,:)), norm(X(j,:) - scrX(x0_index,:)));
            end
        end
    end

    % Construct offset for RHS vector
%     b_offset = zeros(2*K*(K-1),1);
%     for i = 1:K
%         for j = 1:K
%             if i ~= j            
%                 if j < i
%                     row = (i-1)*(K-1) + j;
%                 elseif j > i
%                     row = (i-1)*(K-1) + j - 1;
%                 end
%                 b_offset(row) = ((sqrt(sigma2hat_by_macro(m, i)) + sqrt(sigma2hat_by_macro(m, j)))/sqrt(num_sim))*t_crit;
%                 b_offset(row + K*(K-1)) = b_offset(row); % = ((sqrt(sigma2hat_by_macro(m, i)) + sqrt(sigma2hat_by_macro(m, j)))/sqrt(num_sim))*t_crit;
%             end
%         end
%     end
    
    % Calculate max c value (for fixed alpha)
    tight_c_terms = (A_muhat - A_sigma*t_crit)./(b_RHS/c);
    c_screen(x0_index) = max(tight_c_terms);
    
end

%% Plotting

figure
stairs(sort(c_screen), 1:length(scrX), 'b', 'LineWidth', 2); %(1:length(scrX))/length(scrX)
set(gca, 'FontSize', 14, 'LineWidth', 2)
xlabel('Lipschitz constant $c$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
ylabel('Number of Solutions in $\mathcal{S}$', 'Interpreter', 'latex', 'FontSize', 16, 'FontWeight', 'bold')
title('Sensitivity of Subset Size to Lipschitz Constant', 'FontSize', 16, 'FontWeight', 'bold')

%%

%C = c_screen;

scatter(scrX(:,1), scrX(:,2), [], c_screen, 'filled');
% for k = 1:size(scrX,1)
%     
% rectangle('Position',[scrX(k,1)-C(k)/2  scrX(k,2)-C(k)/2 C(k) C(k)],'facecolor',[0.2 0.8 0.8],'linestyle','none','Curvature',1)
% 
% end
axis square
hold on

plot(X(:,1),X(:,2),'kd','markerfacecolor','k')
plot(scrX(istar,1),scrX(istar,2),'r*')
plot(scrX(istar,1),scrX(istar,2),'ro')

vch = convhull(scrX(:,1),scrX(:,2));
plot(scrX(vch,1),scrX(vch,2),'k-')
xlim([min(scrX(:,1))-5,max(scrX(:,1))+5])
ylim([min(scrX(:,2))-5,max(scrX(:,2))+5])
xlabel('Reorder quantity','interpreter','latex')
ylabel('Order-to quantity','interpreter','latex')