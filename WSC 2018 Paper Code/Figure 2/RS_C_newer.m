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


QVV = sum(frnd(ones([length(n),100000]),repmat(n',1,100000)-1),1);
LV = quantile(QVV,1-alpha);
LVV = quantile(QVV,1-9*alpha/10);
dLVV = LV-LVV;



options = optimoptions(@quadprog,'Algorithm','interior-point-convex','display','off','FunctionTolerance',1/10*abs(dLVV));
[xh0,L,~,~,lvdp] =quadprog(A,a,B,b,[],[],[],[],[],options);
size(lvdp.ineqlin)
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
Ao(K+1:end,K+1:end) = 10^(-12)*eye(length(ao)-K);


%Duality
%2*L+c+lambba*(Bo-bo)

%Ao\(Bo'*lambda)

%Bo*x = bo

%Ao*x+ao+nu*Bo = 0
%nu = Bo\(Ao*x)


xh = zeros(size(ao));
VV = sparse(length(bo),length(Xtry),0);
XVV = zeros(length(ao),length(Xtry));
lamho = [];
for k = 1:length(Xtry)
    
    x0 = Xtry(k,:);
    for dim = 1:d
        Bo(((size(B,1)+K+1):end),dim*K+(1:K)) = -(diag(x0(dim)-Xu(:,dim)));
    end
    options = optimoptions(@quadprog,'Algorithm','interior-point-convex','display','off','FunctionTolerance',1/1000*abs(dLVV));
    
    LB = -inf;
    UB = inf;
    if k > 1
        [xc1 lam1] = EQP(Ao,ao,Bo(AS,:)',bo(AS),0*xc1);
        ASp = AS;
        while min(lam1) < -10^(-8)
            ASp = setdiff(ASp,ASp(lam1 < -10^(-8)));
            [xc1 lam1] = EQP(Ao,ao,Bo(ASp,:)',bo(ASp),0*xc1);
        end
        lam2 = zeros(size(lgp.ineqlin));
        lam2(ASp) = lam1;
        
        xhbp = Ao\(-Bo'*lam2-ao);
        LB = 2*(1/2*xhbp'*Ao*xhbp+ao'*xhbp+lam2'*(Bo*xhbp-bo))+Offset;
        AS = ASp;
        if LB <= LV
            [xc3 lam3] = EQP(Ao,ao,Bo(AS,:)',bo(AS),0*xc1);
            ASp = AS;
            [MVT,W] = max(Bo*xc3-bo);
            while MVT > 10^(-8)
                ASp = union(ASp,W);
                [xc3 lam3] = EQP(Ao,ao,Bo(ASp,:)',bo(ASp),0*xc1);
                [MVT,W] = max(Bo*xc3-bo);
            end
            lam4 = zeros(size(lgp.ineqlin));
            lam4(ASp) = lam3;
            xhbp = Ao\(-Bo'*lam2-ao);
            if max(Bo*xc3-bo) < 10^(-6)
                UB= 2*(1/2*xc3'*Ao*xc3+ao'*xc3)+Offset;
            end
        end
    end
    if LB > LV
        S(k) = LB;
    elseif UB < LV
        S(k) = UB;
    else
        [xc1,L,~,~,lgp]  =quadprog(Ao,ao,Bo,bo,[],[],[],[],0*xh,options);
        AS = find(lgp.ineqlin >10^(-5));
        S(k) = 2*(1/2*xc1'*Ao*xc1+ao'*xc1)+Offset;
    end
    
end
L = (S < LV );


end
