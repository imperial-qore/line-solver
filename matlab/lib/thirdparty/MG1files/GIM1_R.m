function R=GIM1_R(A,Dual,Algor,varargin)
%GIM1_R determines R matrix of a GI/M/1-Type Markov Chain  
%
%   R=GIM1_R(A,Dual,Algor) computes the minimal nonnegative solution to the 
%   matrix equation R = A0 + R A1 + R^2 A2 + R^3 A3 + ... + R^max A_max, 
%   where A = [A0 A1 A2 A3 ... A_max] has m rows and m*max columns and is
%   a nonnegative matrix, with (A0+A1+A2+...+A_max) irreducible and 
%   stochastic.
%
%   The parameter Dual specifies whether the Ramaswami or Bright
%   duals is used. This parameter can take on 3 values:
%
%       'B': Uses the Bright dual
%       'R': Uses the Ramaswami dual
%       'A': Automatic, uses the Bright dual for a positive
%            recurrent chain and the Ramaswami dual otherwise
%
%   The 'A' typically results in the lowest computation time (when CR is
%   used)
%
%   The parameter Algor determines the class of algorithms used to 
%   compute G. This is realized by computing the G matrix of the 
%   Bright/Ramaswami Dual of the Markov chain.
%
%   Four classes are supported: 
%
%       'FI' : Functional Iterations [Neuts]
%       'CR' : Cyclic Reduction [Bini, Meini]
%       'NI' : Newton Iteration [Perez, Telek, Van Houdt]
%       'RR' : Ramaswami Reduction [Bini, Meini, Ramaswami]
%       'IS' : Invariant Subspace [Akar, Sohraby]
%
%   Optional Parameters:
%   
%       See MG1_FI, MG1_CR, MG1_NI, MG1_RR and MG1_IS for the supported
%       optional parameters

m=size(A,1);
dega=size(A,2)/m-1;

% compute invariant vector of A and the drift
% drift > 1: positive recurrent GIM1, drift < 1: transient GIM1
sumA=A(:,dega*m+1:end);
beta=sum(sumA,2);
% beta = (A_maxd)e + (A_maxd + A_maxd-1)e + ... + (Amaxd+...+A1)e
for i=dega-1:-1:1
    sumA=sumA+A(:,i*m+1:(i+1)*m);
    beta=beta+sum(sumA,2);
end
sumA=sumA+A(:,1:m);
theta=stat(sumA);
drift=theta*beta;

if (strcmp(Dual,'R')|(strcmp(Dual,'A') & drift <=1 )) % RAM dual
    % compute the RAM Dual process
    for i=0:dega
        A(:,i*m+1:(i+1)*m)=diag(theta.^(-1))*A(:,i*m+1:(i+1)*m)'*diag(theta);
    end
else % Bright dual
    if (drift > 1) % A -> positive recurrent GIM1
        % compute the Caudal characteristic of A
        eta=GIM1_Caudal(A);
    else % A -> transient GIM1 (=recurrent MG1)
        eta=MG1_Decay(A);
    end    
    % compute invariant vector of A0+A1*eta+A2*eta^2+...+Amax*eta^max
    sumAeta=eta^dega*A(:,dega*m+1:end);
    for i=dega-1:-1:0
        sumAeta=sumAeta+eta^i*A(:,i*m+1:(i+1)*m);
    end
    theta=stat(sumAeta+(1-eta)*eye(m));
    % compute the Bright Dual process
    for i=0:dega
        A(:,i*m+1:(i+1)*m)=eta^(i-1)*diag(theta.^(-1))*A(:,i*m+1:(i+1)*m)'*diag(theta);
    end
end

switch Algor
    case 'FI'
        G=MG1_FI(A,varargin{:});
    case 'CR'
        G=MG1_CR(A,varargin{:});
    case 'NI'
        G=MG1_NI(A,varargin{:});
    case 'RR'
        G=MG1_RR(A,varargin{:});
    case 'IS'
        G=MG1_IS(A,varargin{:});
    otherwise
        error('MATLAB:GIM1_R:InvalidAlgor',...
            'Algorithm ''%s'' is not supported',Algor)
end

if (strcmp(Dual,'R')|(strcmp(Dual,'A') & drift <=1 )) % RAM dual
    R=diag(theta.^(-1))*G'*diag(theta);
else % Bright dual
    R=diag(theta.^(-1))*G'*diag(theta)*eta;
end    