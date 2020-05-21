function [L, S] = RS_L_old(Xtry,X,Y,alpha,sca)
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
CCC = diag(1./b)*B*diag(sigma2hat)*B'*diag(1./b);
sca = 2*max(abs(diag(1./b)*B*muhat')+sqrt(diag(CCC))*norminv(1-alpha/10/2/size(B,1)));
A = spalloc(K,K,K);
A(1:K,1:K) = diag(1./sigma2hat);
a = spalloc(K,1,K);
a(1:K) = -(muhat./sigma2hat)';
LV = quantile(sum(frnd(ones([length(n),1000000]),repmat(n',1,1000000)-1),1),1-alpha);
LV

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
Ao(K+1,K+1) = 10^(-10);
ao = spalloc(K+1,1,K);
ao(1:K) = -((muhat)./sigma2hat)';

%textprogressbar('computing solutions: ');
for k = 1:length(Xtry)
    
    %textprogressbar(k/length(Xtry)*100);
    x0 = Xtry(k,:);
    
    bo(size(B,1)+K+(1:K)) = norms(repmat(x0,K,1)-Xu);   
    
    p =~any(Bo,2);
    
  [x,L,~,~,da] =quadprog(Ao,ao,Bo,sca*bo,[],[],[],[],[],options); 
  
  S(k) = 2*L+sum((muhat.^2)./sigma2hat); 
end
RE = sort(S);
%textprogressbar('done');

L = (S < LV );


end


function N = norms(X)
N = sqrt(sum(X.^2,2));
end
