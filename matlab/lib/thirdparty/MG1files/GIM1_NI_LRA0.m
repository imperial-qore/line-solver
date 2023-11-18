function R=GIM1_R(A,A0hat,Gamma,Dual,varargin)
% Newton Iteration for GI/M/1-Type Markov Chains, exploiting low-rank A0
%
%   R=GIM1_R(A,A0hat,Gamma,Dual) computes the minimal nonnegative solution to the 
%   matrix equation R = A0 + R A1 + R^2 A2 + R^3 A3 + ... + R^N A_N, 
%   where A = [A1 A2 A3 ... AN] has m rows and m*N columns and is
%   a nonnegative matrix. A0 is also nonnegative and can be written as 
%   A0 = Gamma*A0hat, where Gamma and A0hat are mxr and rxm matrices,
%   respectively. (A0+A1+A2+...+AN) is irreducible and stochastic.
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
%       See MG1_NI_LRA0 for the supported optional parameters

A=[Gamma*A0hat A];
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
    A0hat=diag(theta.^(-1))*A0hat';
    Gamma=Gamma'*diag(theta);
    for i=1:dega
        A(:,i*m+1:(i+1)*m)=diag(theta.^(-1))*A(:,i*m+1:(i+1)*m)'*diag(theta);
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
    A0hat=eta^(-1)*diag(theta.^(-1))*A0hat';
    Gamma=Gamma'*diag(theta);
    for i=1:dega
        A(:,i*m+1:(i+1)*m)=eta^(i-1)*diag(theta.^(-1))*A(:,i*m+1:(i+1)*m)'*diag(theta);
    end
end

G=MG1_NI_LRA0(A(:,m+1:end),A0hat,Gamma,varargin{:});

if (strcmp(Dual,'R')|(strcmp(Dual,'A') & drift <=1 )) % RAM dual
    R=diag(theta.^(-1))*G'*diag(theta);
else % Bright dual
    R=diag(theta.^(-1))*G'*diag(theta)*eta;
end    