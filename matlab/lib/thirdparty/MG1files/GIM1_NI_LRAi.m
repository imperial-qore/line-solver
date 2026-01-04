function R=GIM1_R(A0,Ahat,Gamma,Dual,varargin)
% Newton Iteration for GI/M/1-Type Markov Chains, exploiting low-rank Ai,
% for i=1,...,N. 
%
%   R=GIM1_R(A0,Ahat,Gamma,Dual) computes the minimal nonnegative solution to the 
%   matrix equation R = A0 + R A1 + R^2 A2 + R^3 A3 + ... + R^N A_N, 
%   where A0, A1,...,AN are mxm nonnegative matrices such that
%   (A0+A1+A2+...+AN) is irreducible and stochastic.
%   Ai can be written as Ai = Aihat*Gamma, for i = 1,...,N, with Aihat mxr
%   matrices and Gamma an rxm matrix.
%   Ahat = [A1hat A2hat A3hat ... ANhat] has m rows and r*N columns.
%   
%   The parameter Dual specifies whether the Ramaswami or Bright
%   duals is used. This parameter can take on 3 values:
%
%       'B': Uses the Bright dual
%       'R': Uses the Ramaswami dual
%       'A': Automatic, uses the Bright dual for a positive
%            recurrent chain and the Ramaswami dual otherwise
%
%   Optional Parameters:
%   
%       See MG1_NI_LRAi for the supported optional parameters

m=size(A0,1);
r=size(Gamma,1);
dega=size(Ahat,2)/r;
A=zeros(m,m*(dega+1));
A(:,1:m)=A0;
for i=1:dega
    A(:,i*m+1:(i+1)*m)=Ahat(:,(i-1)*r+1:i*r)*Gamma;
end

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
    A0=diag(theta.^(-1))*A0'*diag(theta);
    Gamma=diag(theta.^(-1))*Gamma';
    Ahatnew=zeros(r,m*dega);
    for i=1:dega
        Ahatnew(:,(i-1)*m+1:i*m)=Ahat(:,(i-1)*r+1:i*r)'*diag(theta);
    end
else % Bright dual
    if (drift > 1) % A -> positive recurrent GIM1
        % compute the Caudal characteristic of A
        eta=GIM1_Caudal(A);
    else % A -> transient GIM1 (=recurrent MG1)
        eta=MG1_Decay(A);
    end    
    % compute invariant vector of A0+A1*eta+A2*eta^2+...+A_N*eta^N
    sumAeta=eta^dega*A(:,dega*m+1:end);
    for i=dega-1:-1:0
        sumAeta=sumAeta+eta^i*A(:,i*m+1:(i+1)*m);
    end
    theta=stat(sumAeta+(1-eta)*eye(m));
    % compute the Bright Dual process
    A0=eta^(-1)*diag(theta.^(-1))*A0'*diag(theta);
    Gamma=diag(theta.^(-1))*Gamma';
    Ahatnew=zeros(r,m*dega);
    for i=1:dega
        Ahatnew(:,(i-1)*m+1:i*m)=eta^(i-1)*Ahat(:,(i-1)*r+1:i*r)'*diag(theta);
    end
end

G=MG1_NI_LRAi(A0,Ahatnew,Gamma,varargin{:});

if (strcmp(Dual,'R')|(strcmp(Dual,'A') & drift <=1 )) % RAM dual
    R=diag(theta.^(-1))*G'*diag(theta);
else % Bright dual
    R=diag(theta.^(-1))*G'*diag(theta)*eta;
end    