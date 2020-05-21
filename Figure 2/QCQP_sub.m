function [xh,lam] = QCQP_sub(A,Atb,Ki,nu,R,lam)
%(A x-b)'(A x-b) s.t. x'Lx leq \nu

p = size(R,2);

S = true;
if S
    xht = [AtA+lam*Ki R;
        R' zeros(p)]\[Atb;zeros(p,1)];
    xh =xht(1:end-p);
    if sum(xh'*Ki*xh)<=nu
        maxv = lam;
        minv= 0;
        mf = 1;
    else
        minv = lam;
        mf = 0;
    end
    
    time_2_break = 0;
    if mf == 0
        while mf == 0
            lam = lam*2;
            xht = [AtA+lam*Ki R;
                R' zeros(p)]\[Atb;zeros(p,1)];
            xh =xht(1:(end-p));
            if sum(xh'*Ki*xh)<=nu
                maxv = lam;
                mf = 1;
            end
            if lam > 10^(12)
                time_2_break = 1;
                break;
            end
            
        end
    end
    if time_2_break == 0
        while abs(sum(xh'*Ki*xh)-nu) > 10^(-4)*nu && (maxv-minv) > 10^(-12)
            lam = (maxv-minv)/2+minv;
            xht = [AtA+lam*Ki R;
                R' zeros(p)]\[Atb;zeros(p,1)];
            xh =xht(1:end-p);
            if sum(xh'*Ki*xh)>nu
                minv = lam;
            else
                maxv = lam;
            end
        end
    end
end
end
