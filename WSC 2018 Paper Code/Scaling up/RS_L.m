function [L, S] = RS_L(Xtry,X,Y,alpha)
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
%[2*(abs(diag(1./b)*B*muhat')) sqrt(diag(CCC))*norminv(1-alpha/10/2/size(B,1))]
sca
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
Ao(K+1,K+1) = 10^(-8);
Ao(K+2,K+2) = 10^(-8);
ao(1:K) = -((muhat)./sigma2hat)';
ao(K+1) = -min(muhat)*10^(-8);
LV = quantile(sum(frnd(ones([length(n),100000]),repmat(n',1,100000)-1),1),1-alpha);
options = optimoptions(@quadprog,'Algorithm','interior-point-convex','display','off');

Offset =  sum((muhat.^2)./sigma2hat);

jvec = zeros(1,K+2)';
jvec(end) = -1;
lambdaO = ones(size(bo));
xO = [muhat';min(muhat);0];
nuO = 1;
xV = zeros(1000,length(ao));
nuV = zeros(1000);
lambdaV = zeros(1000,length(bo));


XINDEX = zeros(1000,d);
NLAM = 0;
for k = 1:size(Xu,1)
    x0 = Xtry(kind(k),:);
    bo(size(bo,1)-((K-1):-1:0)) = sca*norms(repmat(x0,K,1)-Xu)+10^(-4);
    
    if NLAM > 0.5
        [RTS,R] = sort(norms(repmat(x0,NLAM,1)-Xtry(1:NLAM,:)));
    end
    
    LB = -inf;
    UB = inf;
    LLV = inf;
    for lv = 1:min(max(10,d),NLAM)
        l = R(lv);
        nup = nuV(l);
        lambdap = lambdaV(l,:)';
        xp = xV(l,:)';
        
        tolVV = 10^(-4);
        ep = tolVV/max(max(Bo));
        IVS = find((lambdap > ep));
        aoo = nup*ao+jvec;
        
        
        rd = nup*Ao*xp+aoo+Bo(IVS,:)'*lambdap(IVS);
        resid = Bo*xp-bo;
        LLVp =norm(rd) + norm(resid.*(resid<0));
        if LLVp <  LLV
            xO = xp;
            lambdaO = lambdap;
            nuO = nup;
            LLV = LLVp;
        end
        
        if max(abs(rd)) > 10^(-5)
            xpp = -(nup*Ao)\(aoo+Bo(IVS,:)'*lambdap(IVS));
            residpp = Bo*xpp-bo;
            LB = 1/2*nup*xpp'*Ao*xpp+aoo'*xpp+lambdap(IVS)'*residpp(IVS)+nup*(Offset/2-LV/2 );
            if LB > 0
                break;
            end
            LLVp =norm(rd) + norm(residpp.*(residpp<0));
            if LLVp <  LLV 
                xO = xpp;
                lambdaO = lambdap;
                nuO = nup;
                LLV = LLVp;
            end
        else
            LB = 1/2*nup*xp'*Ao*xp+aoo'*xp+lambdap(IVS)'*resid(IVS)+nup*(Offset/2-LV/2 );
            if LB > 0
                break;
            end
            LLVp =norm(rd) + norm(resid.*(resid<0));
            if LLVp <  LLV 
                xO = xp;
                lambdaO = lambdap;
                nuO = nup;
                LLV = LLVp;
            end
        end
        
        
        
        if xp'*Ao*xp/2+ao'*xp+ Offset/2-LV/2 < 0
            if  (max(Bo*xp-bo) < 10^(-5))
                UB = min(UB,jvec'*xp);
            end
            if UB < 0
                break;
            end
            LLVp =norm(rd) + norm(resid.*(resid<0));
            if LLVp <  LLV 
                xO = xp;
                lambdaO = lambdap;
                nuO = nup;
                LLV = LLVp;
            end
        else
            mvvv = [muhat';min(muhat);0];
            dxx = xp-[muhat';min(muhat);0];
            aqf = dxx'*Ao*dxx/2 ;
            bqf = dxx'*Ao*mvvv+ao'*dxx;
            cqf =  -LV/2;
            adj_val = (-bqf+sqrt(bqf^2-4*aqf*cqf))/(2*aqf);
            xppp = [muhat';min(muhat);0]+0.95*adj_val*dxx;
            if  (max(Bo*xppp-bo) < 10^(-5))
                UB = min(UB,jvec'*xppp);
            end
            if UB < 0
                break;
            end
            
            residppp = Bo*xppp-bo;
            rddd = nup*xppp'*Ao*xppp+aoo+Bo(IVS,:)'*lambdap(IVS);
            LLVp =norm(rddd) + norm(residppp.*(residppp<0));
            if LLVp <  LLV 
                xO = xppp;
                lambdaO = lambdap;
                nuO = nup;
                LLV = LLVp;
            end
        end
    end
    if LB > 0
        L(k) = 0;
    elseif UB < 0
        L(k) = 1;
    else
        NLAM = NLAM+1;
        XINDEX(NLAM,:) = x0;
        [~, xV(NLAM,:), nuV(NLAM), lambdaV(NLAM,:)] = LPQC_sub(jvec,Bo,bo,Ao,ao,Offset/2-LV/2,xO,nuO,lambdaO);
    end
end

shufV = randsample(1:size(Xtry,1),size(Xtry,1));
for k = shufV
    x0 = Xtry(k,:);
    bo(size(bo,1)-((K-1):-1:0)) = sca*norms(repmat(x0,K,1)-Xu)+10^(-4);
    if NLAM > 0.5
        [RTS,R] = sort(norms(repmat(x0,NLAM,1)-Xtry(1:NLAM,:)));
    end
    LB = -inf;
    UB = inf;
    LLV = inf;
    for lv = 1:min(max(10,d),NLAM)
        l = R(lv);
        nup = nuV(l);
        lambdap = lambdaV(l,:)';
        xp = xV(l,:)';
        
        tolVV = 10^(-4);
        ep = tolVV/max(max(Bo));
        IVS = find((lambdap > ep));
        aoo = nup*ao+jvec;
        
        rd = nup*Ao*xp+aoo+Bo(IVS,:)'*lambdap(IVS);
        resid = Bo*xp-bo;
        LLVp =norm(rd) + norm(resid.*(resid<0));
        if LLVp <  LLV
            xO = xp;
            lambdaO = lambdap;
            nuO = nup;
            LLV = LLVp;
        end
        
        if max(abs(rd)) > 10^(-5)
            xpp = -(nup*Ao)\(aoo+Bo(IVS,:)'*lambdap(IVS));
            residpp = Bo*xpp-bo;
            LB = 1/2*nup*xpp'*Ao*xpp+aoo'*xpp+lambdap(IVS)'*residpp(IVS)+nup*(Offset/2-LV/2 );
            if LB > 0
                break;
            end
            LLVp =norm(residpp.*(residpp<0));
            if LLVp <  LLV 
                xO = xpp;
                lambdaO = lambdap;
                nuO = nup;
                LLV = LLVp;
            end
        else
            LB = 1/2*nup*xp'*Ao*xp+aoo'*xp+lambdap(IVS)'*resid(IVS)+nup*(Offset/2-LV/2 );
            if LB > 0
                break;
            end
            LLVp =norm(rd) + norm(resid.*(resid<0));
            if LLVp <  LLV 
                xO = xp;
                lambdaO = lambdap;
                nuO = nup;
                LLV = LLVp;
            end
        end
        
        if xp'*Ao*xp/2+ao'*xp+ Offset/2-LV/2 < 0
            if  (max(Bo*xp-bo) < 10^(-5))
                UB = min(UB,jvec'*xp);
            end
            if UB < 0
                break;
            end
            LLVp =norm(rd) + norm(resid.*(resid<0));
            if LLVp <  LLV 
                xO = xp;
                lambdaO = lambdap;
                nuO = nup;
                LLV = LLVp;
            end
        else
            mvvv = [muhat';min(muhat);0];
            dxx = xp-[muhat';min(muhat);0];
            aqf = dxx'*Ao*dxx/2 ;
            bqf = dxx'*Ao*mvvv+ao'*dxx;
            cqf =  -LV/2;
            adj_val = (-bqf+sqrt(bqf^2-4*aqf*cqf))/(2*aqf);
            xppp = [muhat';min(muhat);0]+0.95*adj_val*dxx;
            if  (max(Bo*xppp-bo) < 10^(-5))
                UB = min(UB,jvec'*xppp);
            end
            if UB < 0
                break;
            end
            
            residppp = Bo*xppp-bo;
            rddd = nup*xppp'*Ao*xppp+aoo+Bo(IVS,:)'*lambdap(IVS);
            LLVp =norm(rddd) + norm(residppp.*(residppp<0));
            if LLVp <  LLV 
                xO = xppp;
                lambdaO = lambdap;
                nuO = nup;
                LLV = LLVp;
            end
        end
    end
    if LB > 0
        L(k) = 0;
    elseif UB < 0
        L(k) = 1;
    else
        if NLAM < 1000
            NLAM = NLAM+1;
            XINDEX(NLAM,:) = x0;
            [L(k) , xV(NLAM,:), nuV(NLAM), lambdaV(NLAM,:)] = LPQC_sub(jvec,Bo,bo,Ao,ao,Offset/2-LV/2,xO,nuO,lambdaO);
        else
            L(k) = LPQC_sub(jvec,Bo,bo,Ao,ao,Offset/2-LV/2,xO,nuO,lambdaO);
        end
        
    end
end

end

function N = norms(X)
N = sqrt(sum(X.^2,2));
end
