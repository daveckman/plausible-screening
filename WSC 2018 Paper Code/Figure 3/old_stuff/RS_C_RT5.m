function [L, S, fh] = RS_C_RT5(Xtry,X,Y,alpha)
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
A = 10^(-4)*eye(K+d*K);
A(1:K,1:K) =A(1:K,1:K)+diag(1./sigma2hat);
a = spalloc(K+d*K,1,K);
a(1:K) = -((muhat)./sigma2hat)';

%B = normr(B);
options = optimoptions(@quadprog,'Algorithm','interior-point-convex','display','final');
optionsLP = optimoptions(@linprog,'Algorithm','interior-point','display','off','OptimalityTolerance',10^(-5));
LV = quantile(sum(frnd(ones([length(n),100000]),repmat(n',1,100000)-1),1),1-alpha);
[x,~,~,~,lamQP] =quadprog(A,a,B,b,[],[],[],[],[],options);
Offset = sum((muhat.^2)./sigma2hat);
GG = x'*A*x/2+a'*x+Offset/2-LV/2  ;
if GG > 0
    sdasda()
end

Btn = [B;
    eye(K) zeros(K,size(B,2)-K);
    -eye(K) zeros(K,size(B,2)-K)];
btn = [b; x(1:K)+sqrt(LV*sigma2hat(1:K)'); -(x(1:K)-sqrt(LV*sigma2hat(1:K)'))];
UBv = zeros(size(x));
LBv = zeros(size(x));
for c = 1:length(x)
    ova = zeros(size(x));
    ova(c) =1;
    if c >= 1
        [G,xL,lambdaL] =  LP2_sub(ova,Btn,btn,x,ones(size(Btn,1),1));
    else
        [G,xL,lambdaL] =  LP2_sub(ova,Btn,btn,xL,lambdaL);
    end
    if isnan(G)
        LBv(c) = -inf;
    else
        LBv(c) = G;
    end
    if c >= 1
        [G,xU,lambdaU] =  LP2_sub(-ova,Btn,btn,x,ones(size(Btn,1),1));
    else
        [G,xU,lambdaU] =  LP2_sub(-ova,Btn,btn,xU,lambdaU);
    end
    if isnan(G)
        UBv(c) = inf;
    else
        UBv(c) = -G;
    end
end
Small_v = min(UBv(1:K));
Smaller_v = min(LBv(1:K));

Offset = sum((muhat.^2)./sigma2hat);
GG = x'*A*x/2+a'*x+Offset/2-LV/2  ;
if GG > 0
    sdasda()
end
% GG
% [LBv UBv]


x0 = Xtry(1,:);
S = zeros(length(Xtry),1);

Bo = spalloc(K^2+K,K+d*K+2,2*K^2+d*K^2+10*K);
Bo(1:size(B,1),1:size(B,2)) = B;
Bo((size(B,1)+1):(size(B,1)+K),1:K) = -eye(K);
Bo((size(B,1)+1):(size(B,1)+K),end-1) = 1;
Bo((size(B,1)+1):(size(B,1)+K),end) = 1;
Bo((size(B,1)+K+1):(size(B,1)+2*K),1:K) = eye(K);
Bo((size(B,1)+K+1):(size(B,1)+2*K),end-1) = -1;

%-m_0 \leq -m_k - Q
%m_0 \geq m_k +Q
%If Q \geq 0, problem

bo = spalloc(size(Bo,1),1,0);


for dim = 1:d
    Bo(((end-K+1):end),dim*K+(1:K)) = -(diag(x0(dim)-Xu(:,dim)));
end
%Bo = normr(Bo);

p = find(any(Bo,2));
Bo = Bo(p,:);
bo=bo(p)+10^(-4);



Ao = spalloc(K+dim*K+2,K+dim*K+2,K+dim*K+2);
Ao = 10^(-8)*eye(K+dim*K+2);
Ao(1:K,1:K) =Ao(1:K,1:K)+diag(1./sigma2hat);
ao = spalloc(K+dim*K+2,1,K);

ao(1:K) = -((muhat)./sigma2hat)';


jvec = zeros(1,length(ao))';
jvec(end) = 1;

options = optimoptions(@quadprog,'Algorithm','interior-point-convex','display','off');
[xO,~,~,~,lamQP] =  quadprog(A,a,B,b,[],[],[],[],[],options);
xOS = [xO;min(xO(1:K));0];
jvdada = ones(length(bo)-length(b),1);
jvdada(Bo(size(B,1)+1:end,:)*xOS-bo(size(B,1)+1:end) > -10^(-2)) = 10;
jvdada(Bo(size(B,1)+1:end,:)*xOS-bo(size(B,1)+1:end)  < -10^(-2)) = 10^(-1);

nuOS = 1;
lambdaOS = [max(min(lamQP.ineqlin,1000),10^(-3));jvdada];

xO = xOS;
nuO = nuOS;
lambdaO = lambdaOS;

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
A(1:K,1:K) =A(1:K,1:K)+diag(1./sigma2hat);
a = spalloc(K+d*K,1,K);
a(1:K) = -((muhat)./sigma2hat)';

%B = normr(B);
options = optimoptions(@quadprog,'Algorithm','interior-point-convex','display','final');
optionsLP = optimoptions(@linprog,'display','off','OptimalityTolerance',10^(-1),'MaxIterations',100,'ConstraintTolerance',10^(-3));
x =quadprog(A,a,B,b,[],[],[],[],[],options);


LV = quantile(sum(frnd(ones([length(n),1000000]),repmat(n',1,1000000)-1),1),1-alpha);


Offset = sum((muhat.^2)./sigma2hat);
GG = x'*A*x/2+a'*x+Offset/2-LV/2  ;
if GG > 0
    sdasda()
end
GG


x0 = Xtry(1,:);
S = zeros(length(Xtry),1);

Bo = spalloc(K^2+K,K+d*K+2,2*K^2+d*K^2+10*K);
Bo(1:size(B,1),1:size(B,2)) = B;
Bo((size(B,1)+1):(size(B,1)+K),1:K) = -eye(K);
Bo((size(B,1)+1):(size(B,1)+K),end-1) = 1;
Bo((size(B,1)+1):(size(B,1)+K),end) = 1;
Bo((size(B,1)+K+1):(size(B,1)+2*K),1:K) = eye(K);
Bo((size(B,1)+K+1):(size(B,1)+2*K),end-1) = -1;

bo = spalloc(size(Bo,1),1,0);


for dim = 1:d
    Bo(((end-K+1):end),dim*K+(1:K)) = -(diag(x0(dim)-Xu(:,dim)));
end
%Bo = normr(Bo);

p = find(any(Bo,2));
Bo = Bo(p,:);
bo=bo(p)+10^(-4);



Ao = spalloc(K+dim*K+2,K+dim*K+2,K+dim*K+2);
Ao = 10^(-8)*eye(K+dim*K+2);
Ao(1:K,1:K) =Ao(1:K,1:K)+diag(1./sigma2hat);
ao = spalloc(K+dim*K+2,1,K);

ao(1:K) = -((muhat)./sigma2hat)';


jvec = zeros(1,length(ao))';
jvec(end) = -1;

options = optimoptions(@quadprog,'Algorithm','interior-point-convex','display','off');
[xO,~,~,~,lamQP] =  quadprog(A,a,B,b,[],[],[],[],[],options);
xOS = [xO;min(xO(1:K));0];
jvdada = ones(length(bo)-length(b),1);
jvdada(Bo(size(B,1)+1:end,:)*xOS-bo(size(B,1)+1:end) > -10^(-2)) = 10;
jvdada(Bo(size(B,1)+1:end,:)*xOS-bo(size(B,1)+1:end)  < -10^(-2)) = 10^(-1);

nuOS = 1;
lambdaOS = [max(min(lamQP.ineqlin,1000),10^(-3));jvdada];


xO = xOS;
nuO = nuOS;
lambdaO = lambdaOS;

XINDEX = zeros(1000,d);
NLAM = 0;
nel = 0;
shufV = randsample(1:size(Xtry,1),size(Xtry,1));

lambdaV = zeros(1000,2*K);


for k = shufV
    LLV = inf;
    UB = inf;
    LB = -inf;
    x0 = Xtry(k,:);
    for dim = 1:d
        Bo(((end-K+1):end),dim*K+(1:K)) = -(diag(x0(dim)-Xu(:,dim)));
    end
    
    B2 = Bo(((end-K+1):end),1:(end-2));
    
    Q1 = B2.*repmat(LBv',[size(B2,1),1]);
    Q1(isnan(Q1)) =0;
    Q2 =B2.*repmat(UBv',[size(B2,1),1]);
    Q2(isnan(Q2)) =0;
    TT = max((sum(min(Q1,Q2)'))');
    
    xO = xO([1:K end-1 end]);
    LVV = inf;
    if TT < Small_v
        %LB = TT-Small_v;
        L(k) = 1;
        
        B3 = Bo(((end-K+1):end), [K+1:end-2]);
        Q3 = B3.*repmat(LBv(K+1:end)',[size(B3,1),1]);
        Q3(isnan(Q3)) =0;
        Q4 = B3.*repmat(UBv(K+1:end)',[size(B3,1),1]);
        Q4(isnan(Q4)) =0;
        RHS = sum(max(-Q3,-Q4)');
        RHSIDS = find(~isinf(RHS));
        
        Boo = sparse( [Bo(((end-2*K+1):(end-K)), [1:K end-1 end]);
            Bo((end-K+RHSIDS), [1:K end-1 end])]);
        boo = [zeros(K,1);
            RHS(RHSIDS)']+10^(-4);
        
        if NLAM > 0.5
            [RTS,R] = sort(norms(repmat(x0,NLAM,1)-Xtry(1:NLAM,:)));
        end
        
        LLV = inf;
        for lv = 1:min(max(5,d),NLAM)
            l = R(lv);
            nup = nuV(l);
            lambdap = lambdaV(l,[1:K (K+RHSIDS)])';
            xp = xV(l,:)';
            
            resid = Boo*xp-boo;
            aoo = nup*ao([1:K end-1 end])+jvec([1:K end-1 end]);
            Aoo = Ao([1:K end-1 end],[1:K end-1 end]);
            rd = nup*Aoo*xp+aoo+Boo'*lambdap;
            LLVp = norm(rd) + norm(resid.*(resid<0));
            if LLVp <  LLV
                xO = xp;
                lambdaO = lambdap;
                nuO = nup;
                LLV = LLVp;
            end
        end
        
        
        if NLAM < 1000
            if ~isinf(LLV)
                NLAM = NLAM+1;
                XINDEX(NLAM,:) = x0;
                
                [L(k) , xV(NLAM,:), nuV(NLAM), lambdaXX] = LPQC_sub(jvec([1:K end-1 end]),Boo,boo,Ao([1:K end-1 end],[1:K end-1 end]),ao([1:K end-1 end]),Offset/2-LV/2,xO,nuO,lambdaO);
                lambdaV(NLAM,[1:K (K+RHSIDS)])  = lambdaXX;
            else
                NLAM = NLAM+1;
                XINDEX(NLAM,:) = x0;
                [L(k) , xV(NLAM,:), nuV(NLAM), lambdaXX] = LPQC_sub(jvec([1:K end-1 end]),Boo,boo,Ao([1:K end-1 end],[1:K end-1 end]),ao([1:K end-1 end]),Offset/2-LV/2,xOS([1:K end-1 end]),1,ones(size(boo)));
                
                lambdaV(NLAM,[1:K (K+RHSIDS)]) = lambdaXX;
            end
        else
            L(k) = LPQC_sub(jvec([1:K end-1 end]),Boo,boo,Aoo,ao([1:K end-1 end]),Offset/2-LV/2,xO,nuO,lambdaO);
        end
    else
        L(k) = 0;
    end
end
nel


end




function N = norms(X)
N = sqrt(sum(X.^2,2));
end