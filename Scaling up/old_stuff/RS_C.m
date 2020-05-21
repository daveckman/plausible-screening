function [L, S, fh] = RS_C(Xtry,X,Y,alpha)
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
A(1:K,1:K) = diag(1./sigma2hat);
a = spalloc(K+d*K,1,K);
a(1:K) = -((muhat)./sigma2hat)';

Offset = sum((muhat.^2)./sigma2hat);
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

LV = quantile(sum(frnd(ones([length(n),1000/alpha]),repmat(n',1,1000/alpha)-1),1),1-alpha);
xvvv = ones(length(kind),size(Bo,2)+2*size(Bo,1));

kind = kind(randperm(length(kind)));
L = 0.5*ones(size(Xtry,1),1)';
S = 0*ones(size(Xtry,1),1)';


options = optimoptions(@quadprog,'Algorithm','interior-point-convex','display','off','OptimalityTolerance',max(LV/10000,10^(-8)));

for kval = 1:length(kind)
    x0 = Xtry(kind(kval),:);
    
    for dim = 1:d
        Bo(((size(B,1)+K+1):end),dim*K+(1:K)) = -(diag(x0(dim)-Xu(:,dim)));
    end
    if kval == 1
        [setupinfo.Aeq,setupinfo.Beq,setupinfo.lb,setupinfo.ub,setupinfo.flags,setupinfo.options,setupinfo.defaultopt,setupinfo.optionFeedback]=...
            quadprog_setup(Ao,ao,Bo,bo,[],[],[],[],[],options);
    end
    [Lv,xvvv(kval,:)] = QP2_sub(Ao,ao,Bo,bo,setupinfo);
    
    
    S(kind(kval)) = 2*Lv  + Offset;
    L(kind(kval)) = (S(kind(kval)) < LV );
end
J = 1:size(Xtry,1);

dt = delaunayn(Xtry(kind,:));

[mivsS,ivsS] = sort(S(kind));

knind = J(~ismember(J,kind));
knind = knind(randperm(length(knind)));

xss = xvvv(end,:)';
for kval = 1:length(knind)
    x0 = Xtry(knind(kval),:);
    
    
    [t,P] = tsearchn(Xtry(kind,:),dt,x0);
    if ~isnan(t)
        if min(L(kind(dt(t,:)))) > 0.5
            L(knind(kval)) =1;
            S(knind(kval)) = sum(P.*S(kind(dt(t,:))));
        else
            L(knind(kval)) = 0.5;
        end
    end
    
    if L(knind(kval)) < 0.75
        for dim = 1:d
            Bo(((size(B,1)+K+1):end),dim*K+(1:K)) = -(diag(x0(dim)-Xu(:,dim)));
        end
        maxv = min(10,size(Xu,1));
        if ~isnan(t)
            smallerset = unique([ivsS(1:maxv) dt(t,:)]);
        else
            [~,ivsSS] = sort(sum(((repmat(x0,length(kind),1)-Xtry(kind,:))).^2,2));
            smallerset = unique([ivsS(1:maxv) ivsSS(1:min(10+d,size(Xu,1)))']);
        end
        
        gs = length(smallerset);
        
        indkeeps = [];
        for k = 1:(d+1)
            indkeeps = [indkeeps;smallerset'+(k-1)*K];
        end
        indkeeps = [indkeeps;(d+1)*K+1];
        Aop = sparse(Ao(indkeeps,indkeeps));
        aop = sparse(ao(indkeeps));
        Bop = sparse(Bo(:,indkeeps));
        indVV = find((abs(sum(Bop(:,[1:(gs) size(Bop,2)]),2))<10^(-4)).*((sum(abs(Bop(:,[1:(gs) size(Bop,2)])),2))>10^(-4)));
        Bop = sparse(Bop(indVV,:));
        bop = sparse(bo(indVV));
        Offset2 = sum((muhat(smallerset).^2)./sigma2hat(smallerset));
        
        Lvaln = inf;
        
        for l = smallerset
            lval = xvvv(l,(size(Bo,2)+1):(size(Bo,2)+size(Bo,1)))';
            lval = lval(indVV);
            xval = xvvv(l,1:size(Bo,2))';
            xval = xval(indkeeps);
            sval = -((Bop)*xval-bop);
            sval(sval<10^(-3)) = 10^(-3);
            lval(lval<10^(-3)) = 10^(-3);
            
            rd = Aop*xval+aop+(Bop')*lval;
            rp = sval+(Bop)*xval-bop;
            
            Lvalp =norm(rd)+norm(rp);
            if Lvalp<Lvaln
                start_ind = [xval;lval;sval];
            end
        end
        
        Lv2 = QP2_sub(Aop,aop,Bop,bop,setupinfo,start_ind,-Offset2/2+LV/2);
        
        if (2*Lv2+Offset2) > LV
            L(knind(kval)) = 0;
            S(knind(kval)) = (2*Lv2+Offset2);
        end
        
        if L(knind(kval)) > 0.25            
            Lvaln = inf;
            for l = smallerset
                lval = xvvv(l,(size(Bo,2)+1):(size(Bo,2)+size(Bo,1)))';
                xval = xvvv(l,1:size(Bo,2))';                
                sval = -((Bo)*xval-bo);
                sval(sval<10^(-3)) = 10^(-3);
                lval(lval<10^(-3)) = 10^(-3);
                
                rd = Ao*xval+ao+(Bo')*lval;
                rp = sval+(Bo)*xval-bo;
                
                Lvalp =norm(rd)+norm(rp);
                
                if Lvalp<Lvaln
                    start_ind = [xval;lval;sval];
                end
            end
            Lv = QP2_sub(Ao,ao,Bo,bo,setupinfo,start_ind,-Offset/2+LV/2);
            %end
            
            
            S(knind(kval)) = 2*Lv  + Offset;
            L(knind(kval)) = (S(knind(kval)) < LV );
        end
        
    end
    
    
end
