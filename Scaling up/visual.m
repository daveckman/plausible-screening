clear all,close all, clc
d = 2;
K = 500;

dim_order = [5     3     2     9    10     6     8    12     7    11    13     1     4];
scrX = rand(100000,d);
scrX = 0.5 + 4.5*scrX;
dims = dim_order(1:d);

warning off
[IDX, C] = kmeans(scrX,K);
X = C;

num_sim(1:K) = 100;
for k = 1:size(X,1)
    Xsamp = ones(13,1);
    Xsamp(dims) = X(k,:);
    for l = 1:num_sim(k)
        Y(k,l) = SAN(Xsamp,100,round(rand*10^8));
    end
    for l = (num_sim(k)+1):(max(num_sim))
        Y(k,l) = nan;
    end
end
close all

R = RS_C_A(scrX,X,Y,0.95);
plot(scrX(:,1),scrX(:,2),'b.','markersize',4)
hold on
plot(scrX(find(R>0),1),scrX(find(R>0),2),'r.','markersize',2,'markerfacecolor','r')