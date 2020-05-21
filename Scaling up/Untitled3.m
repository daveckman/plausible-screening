close all, clear all, clc
scrX = rand(100000,13);
scrX = 0.5 + 4.5*scrX;


for k = 1:size(scrX,1)
    k
    Y(k) = SAN(scrX(k,:),1000,round(rand*10^8));
    
end

Yc = Y - mean(Y);


Xr =[ scrX scrX.^2];
betas = (Xr'*Xr)^(-1)*(Xr'*Yc');

sigma2 = mean((Yc'-Xr*betas).^2);
sqrt(sigma2)
[betas - 2*sqrt(sigma2*diag((Xr'*Xr)^(-1))) betas + 2*sqrt(sigma2*diag((Xr'*Xr)^(-1)))] 

ME = abs(betas)./sqrt(sigma2*diag((Xr'*Xr)^(-1)));

[~,r]= sort(sum([ME(1:13) ME(14:end)]'));

r

