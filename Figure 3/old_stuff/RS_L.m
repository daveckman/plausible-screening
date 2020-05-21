function [L, S] = RS_L(Xtry,X,Y,alpha,sca)
Xu = unique(X,'rows');
for k =1:length(Xu)
    Ym = Y(ismember(X,Xu(k,:),'rows'));
    n(k) = length(Ym);
    sigma2hat(k) = var(Ym);
    muhat(k) = mean(Ym);
    ES = ismember(Xtry,Xu(k,:),'rows');
    kind(k) = find(ES,1);
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
CCC = diag(1./b)*B*diag(sigma2hat)*B'*diag(1./b);
sca = 2*max(abs(diag(1./b)*B*muhat')+sqrt(diag(CCC))*norminv(1-alpha/10/2/size(B,1)));
sca
A = spalloc(K,K,K);
A(1:K,1:K) = diag(1./sigma2hat);
a = spalloc(K,1,K);
a(1:K) = -(muhat./sigma2hat)';

sca

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
    p = find(any(Bo,2));
    Bo = Bo(p,:);
    bo=bo(p);

Ao = spalloc(K+1,K+1,K);
Ao(1:K,1:K) = diag(1./sigma2hat);
Ao(K+1,K+1) = 10^(-5);
ao = spalloc(K+1,1,K);
ao(1:K) = -((muhat)./sigma2hat)';
LV = quantile(sum(frnd(ones([length(n),100000]),repmat(n',1,100000)-1),1),1-alpha);
options = optimoptions(@quadprog,'Algorithm','interior-point-convex','display','off','OptimalityTolerance',LV/100);

Offset =  sum((muhat.^2)./sigma2hat);


J = 1:size(Xtry,1);

knind = J(~ismember(J,kind));
knind = knind(randperm(length(knind)));
for kval = 1:length(kind)
    x0 = Xtry(kind(kval),:);
    if kval == 1
        [setupinfo.Aeq,setupinfo.Beq,setupinfo.lb,setupinfo.ub,setupinfo.flags,setupinfo.options,setupinfo.defaultopt,setupinfo.optionFeedback]=...
            quadprog_setup(Ao,ao,Bo,bo,[],[],[],[],[],options);
    end
    
    bo(size(bo,1)-((K-1):-1:0)) = norms(repmat(x0,K,1)-Xu);
    
    boo = sca*bo + 10^(-3);
    
    [Lv,xvvv(kval,:)] = QP2_sub(Ao,ao,Bo,boo,setupinfo);
    
    S(kind(kval)) = 2*Lv  + Offset;
    L(kind(kval)) = (S(kind(kval)) < LV );
end


for k = 1:length(Xtry)
    x0 = Xtry(k,:);
    boo = bo; 
    boo(size(bo,1)-((K-1):-1:0)) = norms(repmat(x0,K,1)-Xu);
    boo = sca*boo + 10^(-3);
    
    Lvaln = inf;
    for l = 1:length(kind)
        lval = xvvv(l,(size(Bo,2)+1):(size(Bo,2)+size(Bo,1)))';
        xval = xvvv(l,1:size(Bo,2))';
        
        sval = -((Bo)*xval-boo);
        sval(sval<10^(-3)) = 10^(-3);
        lval(lval<10^(-3)) = 10^(-3);
        
        rd = Ao*xval+ao+(Bo')*lval;
        rp = sval+(Bo)*xval-boo;
                
        Lvalp =norm(rd)+norm(rp);
        if Lvalp<Lvaln
            start_ind = [xval;lval;sval];
        end
    end
    
    
     
    
  [L2,x2] =QP2_sub(Ao,ao,Bo,boo,setupinfo,start_ind,LV/2-Offset/2);
    
    S(k) = 2*L2 +Offset;    
end
%textprogressbar('done');


L = (S < LV );


end


function N = norms(X)
N = sqrt(sum(X.^2,2));
end
