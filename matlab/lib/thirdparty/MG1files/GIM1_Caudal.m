function [eta,v]=GIM1_CAUDAL(A,varargin)
%GIM1_CAUDAL Computes the Spectral Radius of R 
%
%   eta=GIM1_CAUDAL(A) computes the dominant eigenvalue of the
%   matrix R, the smallest nonnegative solution to 
%   R= A0 + R A1 + R^2 A2 + ... + R^max Amax, if the GI/M/1
%   type Markov chain characterized by A is recurrent
%
%   [eta,v]=GIM1_CAUDAL(A) computes the dominant eigenvalue of the
%   matrix R, the smallest nonnegative solution to 
%   R= A0 + R A1 + R^2 A2 + ... + R^max Amax, if the GI/M/1
%   type Markov chain characterized by A is recurrent.
%   eta is the unique solution of PF(A0 + A1 z + A2 z^2 + ... 
%   + Amax z^max) = z on (0,1), where PF denotes the Peron-Frobenius 
%   eigenvalue. The right eigenvector v corresponding to the
%   Peron-Frobenius eigenvalue of A(eta) is also returned.  
%
%
%   Optional Parameters:
%   
%       Dual: When set to 1, the dominant eigenvalue of the Ramaswami
%             dual is returned. The input chain must be a transient
%             M/G/1 type Markov chain.

OptionNames=['Dual        '];
OptionTypes=['numeric'];
OptionValues=[];
 
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.Dual=0;

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

m=size(A,1);
dega=size(A,2)/m-1;

if (options.Dual == 1)
    % the dominant eigenvalue of the Ramaswami dual is
    % identical to the dominant eigenvalue of the GI/M/1 process
    % characterized by A0,A1, ...
end

eta_min=0;
eta_max=1;
eta=1/2;
while (eta_max - eta_min > 10^(-15))
    temp=A(:,dega*m+1:end);
    for i=dega-1:-1:0
        temp=temp*eta+A(:,i*m+1:(i+1)*m);
    end    
    new_eta=max(eig(temp));
    if (new_eta > eta)
        eta_min=eta;
    else
        eta_max=eta;
    end
    eta=(eta_min+eta_max)/2;
end

if (nargout > 1)
    [V,D]=eig(temp);
    v=V(:,find(sum(D,1)==max(sum(D,1))));
end