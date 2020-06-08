function [Lval,xhn,AS] = QPAS_sub(A,b,c,Am,bm,xh,AS)
AS = [];
NAS = [1:length(bm)]';
xhn = 0*xh;
Lim = inf;
inset = 0.5;

bm  = bm + 10^(-3);

dhat = xh;


B = (bm(NAS)-Am(NAS,:)*xhn)./(Am(NAS,:)*dhat);
B(B<10^(-8)) = inf;
[alphak,IM] = min(B);
if ~isempty(B(B>=0))
    B(B<0) = inf;
    [alphak,IM] = min(B);
    
    if alphak > 1
        alphak = (1-10^(-6));
    else
        alphak = alphak*(1-10^(-6));
        AS = [AS;NAS(IM)];
        NAS= NAS([1:(IM-1) (IM+1):end]);
    end
else
    alphak = (1-10^(-6));
end
    xhn = xhn+alphak*dhat;
    
for l = 1:1000
    [dhat,lambda] = EQP(A,b,Am(AS,:)',bm(AS),xhn);
    if ~isempty(lambda)
        while min(lambda) < 0
            [~,i] = min(lambda);
            NAS = [NAS; AS(i)];
            AS = AS([1:(i-1) (i+1):end]);
            lambda= lambda([1:(i-1) (i+1):end]);
        end
    end
    [dhat,lambda] = EQP(A,b,Am(AS,:)',bm(AS),xhn);
    
    
    B = (bm(NAS)-Am(NAS,:)*xhn)./(Am(NAS,:)*dhat);
        B(B<0) = inf;
        [alphak,IM] = min(B);
    if ~isempty(B(B>=0))
        B(B<0) = inf;
        [alphak,IM] = min(B);
        
        if alphak > 1
            alphak = (1-10^(-6));
        else
            alphak = alphak*(1-10^(-6));
            AS = [AS;NAS(IM)];
            NAS= NAS([1:(IM-1) (IM+1):end]);
        end
    else
        alphak = (1-10^(-6));
    end
    
    xhn = xhn+alphak*dhat;
    
    if inset > 0.75 || inset < 0.25
        break;
    end
end
Lval = xhn'*A*xhn + 2*b'*xhn + c;
end
