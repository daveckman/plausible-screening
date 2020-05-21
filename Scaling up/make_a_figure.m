clear all, close all, clc

load('results_07_25_2018')

for i  = 2:7
    subplot(1,6,i-1)
    R = cellstr([num2str(K_vals(1)/1000)]);
    for k = 2:length(K_vals)
        k
        R{k} = [num2str(K_vals(k)/1000)];
    end
    
    boxplot(squeeze(percent_rej(i,:,:))',R)
    title(num2str(d_vals(i)))
    ylim([0,100])
    if i == 1
        
ylabel('Percentage rejected / 1000000')
    end
xlabel(sprintf('Percentage experimented \n ( K / 1000000)'))
end