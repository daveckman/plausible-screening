function [L, S, fh] = RS_C_RT2(Xtry,X,Y,alpha)
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

for k = 1:size(Xu,1)
    x0 = Xu(k,:);
    for dim = 1:d
        Bo(((end-K+1):end),dim*K+(1:K)) = -(diag(x0(dim)-Xu(:,dim)));
    end
    jvdada(Bo(size(B,1)+1:end,:)*xOS-bo(size(B,1)+1:end) > -10^(-2)) = 10;
    jvdada(Bo(size(B,1)+1:end,:)*xOS-bo(size(B,1)+1:end)  < -10^(-2)) = 10^(-1);
    lambdaOS = [max(min(lamQP.ineqlin,1000),10^(-3));jvdada];
    L(k) = LPQC_sub(jvec,Bo,bo,Ao,ao,Offset/2-LV/2,xO,nuO,lambdaO);
    
end


Offset =  sum((muhat.^2)./sigma2hat);

XINDEX = zeros(1000,d);
NLAM = 0;

nel = 0;
shufV = randsample(1:size(Xtry,1),size(Xtry,1));
for k = shufV
    x0 = Xtry(k,:);
    for dim = 1:d
        Bo(((end-K+1):end),dim*K+(1:K)) = -(diag(x0(dim)-Xu(:,dim)));
    end
    if NLAM > 0.5
        [RTS,R] = sort(norms(repmat(x0,NLAM,1)-Xtry(1:NLAM,:)));
    end
    LLV = inf;
    UB = inf;
    for lv = 1:min(max(1,d),NLAM)
        l = R(lv);
        nup = nuV(l);
        lambdap = lambdaV(l,:)';
        xp = xV(l,:)';
        resid = Bo*xp-bo;
        
        aoo = nup*ao+jvec;
        
        rd = nup*Ao*xp+aoo+Bo'*lambdap;
        
        if xp'*Ao*xp/2+ao'*xp+ Offset/2-LV/2 < 0
            if  (max(Bo*xp-bo) < 10^(-5))
                UB = min(UB,jvec'*xp);
            end
            LLVp =norm(rd) + norm(resid.*(resid<0));
            if LLVp <  LLV
                xO = xp;
                lambdaO = lambdap;
                nuO = nup;
                LLV = LLVp;
            end
        else
            if GG <0
                mvvv = xOS;
                dxx = xp-mvvv;
                aqf = dxx'*Ao*dxx/2 ;
                bqf = dxx'*Ao*mvvv+ao'*dxx;
                cqf = mvvv'*Ao*mvvv/2+mvvv'*ao+ Offset/2-LV/2;
                
                adj_val = (-bqf+sqrt(bqf^2-4*aqf*cqf))/(2*aqf);
                xppp = mvvv+adj_val*dxx;
                if  (max(Bo*xppp-bo) < 10^(-5))
                    UB = min(UB,jvec'*xppp);
                end
            end
            
            residppp = Bo*xppp-bo;
            rddd = nup*xppp'*Ao*xppp+aoo+Bo'*lambdap;
            LLVp =norm(rddd) + norm(residppp.*(residppp<0));
            if LLVp <  LLV
                xO = xppp;
                lambdaO = lambdap;
                nuO = nup;
                LLV = LLVp;
            end
        end
        if UB < 0
            break;
        end
    end
    
    
    
    %         LB_LB =  1/2*nup*xp'*Ao*xp+aoo'*xp+nup*(Offset/2-LV/2 );
    %
    %         if LB_LB > 0
    %             [lambdap2,LTT] = linprog(-resid,[],[],Bo',-nup*Ao*xp-aoo,zeros(size(resid)),[],optionsLP );
    %             LTT
    %             if ~isempty(lambdap2)
    %                 LB = LB_LB -LTT;
    %                 lambdap = lambdap2;
    %                 LLVp =norm(resid.*(resid<0));
    %                 if LLVp <  LLV
    %                     xO = xp;
    %                     lambdaO = lambdap;
    %                     nuO = nup;
    %                     LLV = LLVp;
    %                 end
    %             end
    %         end
    %
    %         if LB > 0
    %             break;
    %         end
    
    %   end
    %     if LB > 0
    %         L(k) = 0;
    %         nel = nel +1;
    %     [LB UB]
    %     else
    if UB < 0
        L(k) = 1;
        nel = nel +1;
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
sum(max(lambdaV)>10^(-1))
nel

end




function N = norms(X)
N = sqrt(sum(X.^2,2));
end
