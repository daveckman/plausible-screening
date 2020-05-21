close all, clear all, clc
dim_order = [5     3     2     9    10     6     8    12     7    11    13     1     4];
d_vals = [2 3 4 5 6 7 8 13];
K_vals = [100 250 500 1000 2500 5000];
for i = 1:length(d_vals)
    i
    for j = 1:length(K_vals)
        j
        for reps = 1:10
            d = d_vals(i);
            K = K_vals(j);
            
            scrX = rand(100000,d);
            scrX = 0.5 + 4.5*scrX;
            dims = dim_order(1:d);
            
            warning off
            X = scrX(1:K,:);
            
            num_sim(1:K) = 100;
            for k = 1:size(X,1)
                Xsamp = ones(13,1);
                Xsamp(dims) = X(k,:);
                for l = 1:num_sim(k)
                    Y(k,l) = SAN( Xsamp,100,round(rand*10^8));
                end
                for l = (num_sim(k)+1):(max(num_sim))
                    Y(k,l) = nan;
                end
            end
            
            R = RS_C_A(scrX,X,Y,0.95);
            mean(R<0.5)
            percent_rej(i,j,reps) = 100*mean(R<0.5);
            save('results_07_25_2018','percent_rej','d_vals','K_vals','i','j','reps')
        end
    end
end


disp('herewego')
profile on
tic
toc
sum(R)
profile viewer
