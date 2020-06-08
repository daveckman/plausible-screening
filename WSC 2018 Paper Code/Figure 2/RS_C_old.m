function [L, S, fh] = RS_C(Xtry,X,Y,alpha)
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

B = spalloc(K^2,K+d*K,2*K^2+d*K^2);
for k = 1:K
    B(((k-1)*K+1):(k*K),1:K) = -eye(K);
    B(((k-1)*K+1):(k*K),k) = 1;
    B((k-1)*K+k,k) = 0;
    for dim = 1:d
        B(((k-1)*K+1):(k*K),dim*K+k) = (Xu(k,dim)-Xu(:,dim));
        
    end
end
b = spalloc(size(B,1),1,0);


A = spalloc(K+d*K,K+d*K,K);
A(1:K,1:K) = diag(1./sigma2hat);
a = spalloc(K+d*K,1,K);
a(1:K) = -((muhat)./sigma2hat)';

options = optimoptions(@quadprog,'Algorithm','interior-point-convex','display','off','OptimalityTolerance',10^(-8));
[x,L] =quadprog(A,a,B,b,[],[],[],[],[],options);
Offset = sum((muhat.^2)./sigma2hat);
L = 2*L + Offset;

x0 = Xtry(1,:);
S = zeros(length(Xtry),1);

Bo = spalloc(K^2+2*K,K+d*K+1,2*K^2+d*K^2+10*K);
Bo(1:size(B,1),1:size(B,2)) = B;

Bo((size(B,1)+1):(size(B,1)+K),1:K) = -eye(K);
Bo((size(B,1)+1):(size(B,1)+K),end) = 1;
Bo((size(B,1)+K+1):(size(B,1)+2*K),1:K) = eye(K);
Bo((size(B,1)+K+1):(size(B,1)+2*K),end) = -1;

for dim = 1:d
    Bo(((size(B,1)+K+1):end),dim*K+(1:K)) = -(diag(x0(dim)-Xu(:,dim)));
end
bo = spalloc(size(Bo,1),1,0);

Ao = spalloc(K+dim*K+1,K+dim*K+1,K);
Ao(1:K,1:K) = diag(1./sigma2hat);
ao = spalloc(K+dim*K+1,1,K);
ao(1:K) = -((muhat)./sigma2hat)';
Ao(K+1:end,K+1:end) = 10^(-6)*eye(length(ao)-K);
x = zeros(length(ao));
for k = 1:length(Xtry)
    x0 = Xtry(k,:);
    for dim = 1:d
        Bo(((size(B,1)+K+1):end),dim*K+(1:K)) = -(diag(x0(dim)-Xu(:,dim)));
    end
    [x,L] =quadprog(Ao,ao,Bo,bo,[],[],[],[],[],options,x);
    S(k) = 2*L  + Offset;
end
S
LV = quantile(sum(frnd(ones([length(n),100000]),repmat(n',1,100000)-1),1),1-alpha);
L = (S < LV );


end
