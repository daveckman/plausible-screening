function x = lsqnonneg_me(C,d)
%LSQNONNEG Linear least squares with nonnegativity constraints.
%   X = LSQNONNEG(C,d) returns the vector X that minimizes NORM(d-C*X)
%   subject to X >= 0. C and d must be real.
%
%   X = LSQNONNEG(C,d,OPTIONS) minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, an argument
%   created with the OPTIMSET function.  See OPTIMSET for details. Used
%   options are Display and TolX. (A default tolerance TolX of 
%   10*MAX(SIZE(C))*NORM(C,1)*EPS is used.) 
%   
%   X = LSQNONNEG(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with the matrix 'C' in PROBLEM.C, the vector 'd' in
%   PROBLEM.d, the options structure in PROBLEM.options, and solver name
%   'lsqnonneg' in PROBLEM.solver. 
%
%   [X,RESNORM] = LSQNONNEG(...) also returns the value of the squared 2-norm of 
%   the residual: norm(d-C*X)^2.
%
%   [X,RESNORM,RESIDUAL] = LSQNONNEG(...) also returns the value of the  
%   residual: d-C*X.
%   
%   [X,RESNORM,RESIDUAL,EXITFLAG] = LSQNONNEG(...) returns an EXITFLAG that
%   describes the exit condition. Possible values of EXITFLAG and the
%   corresponding exit conditions are
%
%    1  LSQNONNEG converged with a solution X.
%    0  Iteration count was exceeded. Increasing the tolerance
%       (OPTIONS.TolX) may lead to a solution.
%  
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT] = LSQNONNEG(...) returns a structure
%   OUTPUT with the number of steps taken in OUTPUT.iterations, the type of 
%   algorithm used in OUTPUT.algorithm, and the exit message in OUTPUT.message.
%
%   [X,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA] = LSQNONNEG(...) returns 
%   the dual vector LAMBDA  where LAMBDA(i) <= 0 when X(i) is (approximately) 0 
%   and LAMBDA(i) is (approximately) 0 when X(i) > 0.
% 
%   See also LSCOV, SLASH.

%   Copyright 1984-2016 The MathWorks, Inc. 

% Reference:
%  Lawson and Hanson, "Solving Least Squares Problems", Prentice-Hall, 1974.

% Check if more inputs have been passed. In that case error.

n = size(C,2);
% Initialize vector of n zeros and Infs (to be used later)
nZeros = zeros(n,1);
wz = nZeros;

P = false(n,1);
Z = true(n,1);
x = nZeros;

resid = d - C*x;
w = C'*resid;

% Set up iteration criterion
outeriter = 0;
iter = 0;
itmax = floor(n/3);

tol = 10^(-5);
% Outer loop to put variables into set to hold positive coefficients
while any(Z) && any(w(Z) > tol)  && max(abs(resid)) > 10^(-2)
    if itmax < outeriter
        break;
    end
   outeriter = outeriter + 1;
   % Reset intermediate solution z
   z = nZeros; 
   % Create wz, a Lagrange multiplier vector of variables in the zero set.
   % wz must have the same size as w to preserve the correct indices, so
   % set multipliers to -Inf for variables outside of the zero set.
   wz(P) = -Inf;
   wz(Z) = w(Z);
   % Find variable with largest Lagrange multiplier
   [~,t] = max(wz);
   % Move variable t from zero set to positive set
   P(t) = true;
   Z(t) = false;
   % Compute intermediate solution using only variables in positive set
   z(P) = C(:,P)\d;
   % inner loop to remove elements from the positive set which no longer belong
   while any(z(P) <= 0)
       iter = iter + 1;
       % Find indices where intermediate solution z is approximately negative
       Q = (z <= 0) & P;
       % Choose new x subject to keeping new x nonnegative
       alpha = min(x(Q)./(x(Q) - z(Q)));
       x = x + alpha*(z - x);
       % Reset Z and P given intermediate values of x
       Z = ((abs(x) < tol) & P) | Z;
       P = ~Z;
       z = nZeros;           % Reset z
       z(P) = C(:,P)\d;      % Re-solve for z
   end
   x = z;
   resid = d - C*x;
   w = C'*resid;
end
x = x.*(x>=0);
