function [L, S] = RS_L_RT(Xtry,X,Y,alpha,sca)
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

B = spalloc(K^2,K,4*K^2);
b = ones(size(B,1),1);
for k = 1:K
    B(((k-1)*K+1):(k*K),1:K) = -eye(K);
    B(((k-1)*K+1):(k*K),k) = 1;
    B(((k-1)*K+k),k) = 0;
    b(((k-1)*K+1):(k*K)) = norms(repmat(Xu(k,:),K,1)-Xu);
end
CCC = diag(1./b)*B*diag(sigma2hat)*B'*diag(1./b);
sca = 2*max(abs(diag(1./b)*B*muhat')+sqrt(diag(CCC))*norminv(1-alpha/10/2/size(B,1)));
b = sca*b;
A = spalloc(K,K,K);
A(1:K,1:K) = diag(1./sigma2hat);
a = spalloc(K,1,K);
a(1:K) = -(muhat./sigma2hat)';


S = zeros(length(Xtry),1);
Bo = spalloc(2*K^2+2*K,K+2,4*K^2+4*K);
bo = ones(size(Bo,1),1);
bo(1:size(B,1)) = b;
Bo(1:size(B,1),1:size(B,2)) = B;
Bo(size(B,1)+(1:K),1:K) = -eye(K);
Bo(size(B,1)+(1:K),K+1) = 1;
Bo(size(B,1)+(1:K),K+2) = 1;
bo(size(B,1)+(1:K)) = 0;
Bo(size(B,1)+K+(1:K),1:K) = eye(K);
Bo(size(B,1)+K+(1:K),K+1) = -1;
    p = find(any(Bo,2));
    Bo = Bo(p,:);
    bo=bo(p);




Ao = spalloc(K+2,K+2,K);
Ao(1:K,1:K) = diag(1./sigma2hat);
ao = spalloc(K+2,1,K);
Ao(K+1,K+1) = 10^(-4);
Ao(K+2,K+2) = 10^(-4);
ao(1:K) = -((muhat)./sigma2hat)';
LV = quantile(sum(frnd(ones([length(n),100000]),repmat(n',1,100000)-1),1),1-alpha);
options = optimoptions(@quadprog,'Algorithm','interior-point-convex','display','off');




Offset =  sum((muhat.^2)./sigma2hat);

jvec = zeros(1,K+2)';
jvec(end) = -1;

LV = quantile(sum(frnd(ones([length(n),100000]),repmat(n',1,100000)-1),1),1-alpha);


lambdaO = ones(size(bo));
xO = [muhat';min(muhat);0];
nuO = 1;
xV = zeros(size(Xu,1),length(ao));
nuV = zeros(size(Xu,1));
lambdaV = zeros(size(Xu,1),length(bo));
for k = 1:size(Xu,1)
    x0 = Xu(k,:);
    bo(size(bo,1)-((K-1):-1:0)) = sca*norms(repmat(x0,K,1)-Xu)+10^(-4);
    

    [~, xV(k,:), nuV(k), lambdaV(k,:)] = LPQC_sub(jvec,Bo,bo,Ao,ao,Offset/2-LV/2,xO,nuO,lambdaO);
    nuO = nuV(k);
    
    quadV(k,:) = Ao*( xV(k,:)');%xV(k,:)*
    
    
end

tt1 = 0;
tt2 =0;
for k = 1:size(Xtry,1)
    x0 = Xtry(k,:);
    bo(size(bo,1)-((K-1):-1:0)) = sca*norms(repmat(x0,K,1)-Xu)+10^(-4);
    
    [RTS,R] = sort(norms(repmat(x0,K,1)-Xu));
    
    Lvaln = inf;
    trav1 = tic;
    LB = -inf;
    UB = inf;
    for lv = 1:8
        l = R(lv);
        xp= xV(l,:)';        
        nup = nuV(l);
        lambdap = lambdaV(l,:)';
        sp = (-Bo*xp+bo);
        sp(sp <= 10^(-5)) = 10^(-5);
        tolVV = 10^(-4);
        ep = tolVV/max(max(Bo));
        IVS = find((lambdap > ep));
        
        aoo = nup*ao+jvec;

        rd = nup*quadV(l,:)'+aoo+Bo(IVS,:)'*lambdap(IVS);
        rp = sp+Bo*xp-bo;
                
        Lvalp =norm(rd)+norm(rp);
        if Lvalp<Lvaln
            xs = xp;
            nus = nup;
            lambdas = lambdap;
        end
        resid = Bo*xp-bo;
        
        if max(abs(rd)) < 10^(-3)  && max(abs(rp)) < 10^(-3)
            LB = 1/2*nup*xp'*quadV(l,:)'+aoo'*xp+lambdap'*resid+nup*(Offset/2-LV/2 );
            if LB > 0
            break;
            end
        end
        if  (max(Bo*xp-bo) < 10^(-3)) && (xp'*quadV(l,:)'/2+ao'*xp+ Offset/2-LV/2 < 0)
            UB = jvec'*xp;
            if UB < 0
            break;
            end
        end
    end
    
    tt1 = tt1 + toc(trav1);
    
    if UB < 0
        L(k) = 1;
    elseif LB > 0
        L(k) = 0;
    else        
        trav2 = tic;
        size(ao)
        size(jvec)
        L(k)= LPQC_sub(jvec,Bo,bo,Ao,ao,Offset/2-LV/2,xs,nus,lambdas);
        tt2 = tt2 + toc(trav2);
    end
end
end


function N = norms(X)
N = sqrt(sum(X.^2,2));
end
