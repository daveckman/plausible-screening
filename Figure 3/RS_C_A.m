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
for c = (K+1):(K+d*K)
    kvs = mod(c,K);
    dims = (c-kvs)/K;
    if kvs == 0
        kvs = K;
        dims = c/K-1;
    end
    
    ova2 = zeros(1,d);
    ova2(dims) = 1;
    Avva = (repmat(Xu(kvs,:),[K-1,1])-Xu([1:(kvs-1) (kvs+1):K],:));
    bvva = UBv([1:(kvs-1) (kvs+1):K])-LBv(kvs);
        
    G =  LP_sub(ova2',Avva,bvva,zeros(d,1),ones(size(Avva,1),1));
    if isnan(G)
        LBv(c) = -inf;
    else
        LBv(c) = G;
    end    
    G =  LP_sub(-ova2',Avva,bvva,zeros(d,1),ones(size(Avva,1),1));
    if isnan(G)
        UBv(c) = inf;
    else
        UBv(c) = -G;
    end
    
end

Small_v = min(UBv(1:K));

shufV = 1:size(Xtry,1);
for k = shufV
    x0 = Xtry(k,:);
    
    UB_now = LBv(1:K)';
    for dim = 1:d
        UB_now =  UB_now + ...
            min(-(x0(dim)-Xu(:,dim)).*LBv(dim*K+(1:K)),-(x0(dim)-Xu(:,dim)).*UBv(dim*K+(1:K)))';
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