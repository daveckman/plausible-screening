function [L, S] = RS_L(Xtry,X,Y,alpha,sca)
Xu = unique(X,'rows');
for k =1:length(Xu)
    Ym = Y(ismember(X,Xu(k,:),'rows'));
    n(k) = length(Ym);
    sigma2hat(k) = var(Ym);
    muhat(k) = mean(Ym);
end

spos = sigma2hat(sigma2hat>(10^(-10)*max(sigma2hat)));
sigma2hat(sigma2hat<=quantile(spos,0.1)) = quantile(spos,0.1);

sigma2hat = sigma2hat./n;

K = size(Xu,1);
d = size(Xu,2);

B = spalloc(2*K^2,K,4*K^2);
b = ones(size(B,1),1);
for k = 1:K
    B(((k-1)*K+1):(k*K),1:K) = -eye(K);
    B(((k-1)*K+1):(k*K),k) = 1;
    B(((k-1)*K+k),k) = 0;
    b(((k-1)*K+1):(k*K)) = norms(repmat(Xu(k,:),K,1)-Xu);
    B((K^2 + (k-1)*K+1):(K^2 + k*K),1:K) = eye(K);
    B((K^2 + (k-1)*K+1):(K^2 + k*K),k) = -1;
    B(K^2+((k-1)*K+k),k) = 0;
    b(K^2+(((k-1)*K+1):(k*K))) = norms(repmat(Xu(k,:),K,1)-Xu);
end
A = spalloc(K,K,K);
A(1:K,1:K) = diag(1./sigma2hat);
a = spalloc(K,1,K);
a(1:K) = -(muhat./sigma2hat)';

options = optimoptions(@quadprog,'Algorithm','interior-point-convex','display','off');
S = zeros(length(Xtry),1);
Bo = spalloc(2*K^2+2*K,K+1,4*K^2+4*K);
bo = ones(size(Bo,1),1);
bo(1:size(B,1)) = b;
Bo(1:size(B,1),1:size(B,2)) = B;
Bo(size(B,1)+(1:K),1:K) = -eye(K);
Bo(size(B,1)+(1:K),end) = 1;
bo(size(B,1)+(1:K)) = 0;
Bo(size(B,1)+K+(1:K),1:K) = eye(K);
Bo(size(B,1)+K+(1:K),end) = -1;

Ao = spalloc(K+1,K+1,K);
Ao(1:K,1:K) = diag(1./sigma2hat);
ao = spalloc(K+1,1,K);
ao(1:K) = -((muhat)./sigma2hat)';

for k = 1:length(Xtry)
    x0 = Xtry(k,:);
    
    bo(size(B,1)+K+(1:K)) = norms(repmat(x0,K,1)-Xu);    
    [x,L] =quadprog(Ao,ao,Bo,sca*bo,[],[],[],[],[],options);    
    S(k) = 2*L + sum((muhat.^2)./sigma2hat);    
end

LV = quantile(sum(frnd(ones([length(n),100000]),repmat(n',1,100000)-1),1),1-alpha);

L = (S < LV );


end


function N = norms(X)
N = sqrt(sum(X.^2,2));
end
