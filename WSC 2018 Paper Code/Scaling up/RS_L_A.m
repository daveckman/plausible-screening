function [L, S, fh] = RS_C_A(Xtry,X,Y,alpha)
Xu = X;
for k =1:size(Xu,1)
    if sum(isnan(Y(k,:))) > 0
        ne = find(isnan(Y(k,:)),1,'first');
        Ym = Y(k,1:(ne-1));
    else
        Ym = Y(k,:);
    end
    n(k) = length(Ym);
    sigma2hat(k) = var(Ym);
    muhat(k) = mean(Ym);
end

spos = sigma2hat(sigma2hat>(10^(-10)*max(sigma2hat)));
sigma2hat(sigma2hat<=quantile(spos,0.1)) = quantile(spos,0.1);

sigma2hat = sigma2hat./n;

K = size(Xu,1);
d = size(Xu,2);

UBv = zeros(K+d*K,1);
LBv =  zeros(K+d*K,1);
LBv(1:K) = muhat(1:K)'+tinv(alpha/2/K,n'-1).*sqrt(sigma2hat(1:K)');
UBv(1:K) = muhat(1:K)'+tinv(1-alpha/2/K,n'-1).*sqrt(sigma2hat(1:K)');

JUBv = sort(UBv);
Small_v = JUBv(2);

shufV = 1:size(Xtry,1);
for k = shufV
    x0 = Xtry(k,:);
    
    UB_now = LBv(1:K)';
    
    
    UB_now = UB_now-sca*norms(repmat(x0,K,1)-Xu);   
    
    TT = max(UB_now);
    if TT < Small_v
        L(k) = 1;
    else
        L(k) = 0;
    end
end


end




function N = norms(X)
N = sqrt(sum(X.^2,2));
end