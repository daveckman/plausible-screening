function [L, S] = RS_L_CRN(Xtry,X,Y,alpha)
muhat = mean(Y');
Smat = cov(Y');


K = size(X,1);
d = size(X,2);
n = size(Y,2);
Sbar = (Smat+min(diag(Smat))/100*eye(K))/n;


LV = K*(n-1)/(n-K)*finv(1-alpha,n,n-K);
LV

B = spalloc(2*K^2,K,4*K^2);
b = ones(size(B,1),1);
for k = 1:K
    B(((k-1)*K+1):(k*K),1:K) = -eye(K);
    B(((k-1)*K+1):(k*K),k) = 1;
    B(((k-1)*K+k),k) = 0;
    b(((k-1)*K+1):(k*K)) = norms(repmat(X(k,:),K,1)-X);
    B((K^2 + (k-1)*K+1):(K^2 + k*K),1:K) = eye(K);
    B((K^2 + (k-1)*K+1):(K^2 + k*K),k) = -1;
    B(K^2+((k-1)*K+k),k) = 0;
    b(K^2+(((k-1)*K+1):(k*K))) = norms(repmat(X(k,:),K,1)-X);
end
CCC = diag(1./b)*B*Sbar*B'*diag(1./b);
sca = 2*max(abs(diag(1./b)*B*muhat')+sqrt(diag(CCC))*norminv(1-alpha/10/2/size(B,1)));
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
 %   p = find(any(Bo,2));
 %   Bo = Bo(p,:);
  %  bo=bo(p);

Ao = spalloc(K+1,K+1,K);
Ao(1:K,1:K) = Sbar^(-1);
Ao(K+1,K+1) = 10^(-10);
ao = spalloc(K+1,1,K);
ao(1:K) = -Sbar\(muhat');
options = optimoptions(@quadprog,'Algorithm','interior-point-convex','display','off');

Offset =  -muhat*ao(1:K);
for k = 1:length(Xtry)
    x0 = Xtry(k,:);
    
    bo(size(B,1)+K+(1:K)) = norms(repmat(x0,K,1)-X);   
    
    p =~any(Bo,2);
    
  [x,L,~,~,da] =quadprog(Ao,ao,Bo,sca*bo,[],[],[],[],[],options); 
  S(k) = 2*L+Offset; 
end
RE = sort(S);
S(end-10:end)
L = (S < LV );



end


function N = norms(X)
N = sqrt(sum(X.^2,2));
end
