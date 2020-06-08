
function C=oracle_queue(X)
for k = 1:length(X)


c = X(k);
barS = 1000;
mu = 1/barS;
rho = 1*barS;
a = rho/c;
expon = c*mu*(1-a);
expon2 = mu;

B=erlangb(c,rho);
C(k)= X(k)*0.01+(1-erlangc(c,rho))*sqrt(pi)/2/sqrt(expon2)+erlangc(c,rho)*sqrt(pi)/2*(expon2/sqrt(expon)-expon/sqrt(expon2))/(expon2-expon);

end

end


function B=erlangb(n,rho)
  if ((floor(n) ~= n) || (n < 1))
    warning('n is not a positive integer');
    B=NaN;
    return;
  end
  if (rho < 0.0)
    warning('rho is negative!');
    B=NaN;
    return;
  end
B=1;
for k=1:n
  B=((rho*B)/k)/(1+rho*B/k); 
end

end


function C=erlangc(n,rho)
%
% Sanity check- make sure that n is a positive integer.
%
  if ((floor(n) ~= n) || (n < 1))
    warning('n is not a positive integer');
    C=NaN;
    return;
  end;
%
% Sanity check- make sure that rho >= 0.0.
%
  if (rho < 0.0)
    warning('rho is negative!');
    C=NaN;
    return;
  end;
B=erlangb(n,rho);
C=n*B/(n-rho*(1-B));
end


