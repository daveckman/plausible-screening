
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