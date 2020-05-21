function [L, S, fh] = RS_L2_A(Xtry,X,Y,alpha,sca)
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

LBv(1:K) = muhat(1:K)'+tinv(alpha/2/K,n'-1).*sqrt(sigma2hat(1:K)');
UBv(1:K) = muhat(1:K)'+tinv(1-alpha/2/K,n'-1).*sqrt(sigma2hat(1:K)');
sca= 0;
for k = 1:K
    x0 = Xu(k,:);
    for k2 = 1:(k-1)
        x1 = Xu(k2,:);
        abdist = norms(x0-x1);
        INDS = [1:(k2-1) (k2+1):(k-1) (k+1):K];
        acdist = norms(repmat(x0,K-2,1)-Xu(INDS,:))';
        bcdist = norms(repmat(x1,K-2,1)-Xu(INDS,:))';
        max_diff = abs(LBv(k)-(UBv(k2)./(1+abdist./acdist)+UBv(INDS)./(1+acdist./abdist)));
       
        scap = max(max_diff./(bcdist./(1/abdist+1./acdist)));
        if scap > sca
            sca = scap;
        end
    end
    
end

K = size(Xu,1);
d = size(Xu,2);

JUBv = sort(UBv);
Small_v = JUBv(2);

sca
shufV = 1:size(Xtry,1);
for k = shufV
    x0 = Xtry(k,:);
    UB_now = zeros(K,1);
    
   
    for k2 = 1:K
        x1 = Xu(k2,:);
        abdist = norms(x0-x1);
        INDS = [1:(k2-1) (k2+1):K];
        acdist = norms(repmat(x0,K-1,1)-Xu(INDS,:))';
        bcdist = norms(repmat(x1,K-1,1)-Xu(INDS,:))';
        UB_now(k2) = max(-sca*bcdist./(1/abdist+1./acdist)+LBv(k2)./(1+abdist./acdist)+LBv(INDS)./(1+acdist./abdist));
    end
    
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