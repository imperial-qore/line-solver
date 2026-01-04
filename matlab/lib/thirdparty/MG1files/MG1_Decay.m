function [eta,uT]=MG1_Decay(A)
%MG1_Decay Computes the Decay Rate of the MG1 type MC [Falkenberg]
%
%   eta=MG1_Decay(A) computes the decay rate of a recurrent M/G/1 
%   type Markov chain, it is the unique solution of 
%   PF(A0 + A1 z + A2 z^2 + ... + Amax z^max) = z on (1,RA), where 
%   PF denotes the Peron-Frobenius eigenvalue. 
%
%   [eta,uT]=MG1_Decay(A) computes the decay rate of a recurrent M/G/1 
%   type Markov chain and the left eigenvector corresponding to the
%   Peron-Frobenius eigenvalue of A(eta).

m=size(A,1);
dega=size(A,2)/m-1;

eta=1;
new_eta=0;
while (new_eta - eta < 0)
    eta=eta+1;
    temp=A(:,dega*m+1:end);
    for i=dega-1:-1:0
        temp=temp*eta+A(:,i*m+1:(i+1)*m);
    end    
    new_eta=max(eig(temp));
end

eta_min=eta-1;
eta_max=eta;
eta=eta_min+1/2;
while (eta_max - eta_min > 10^(-15))
    temp=A(:,dega*m+1:end);
    for i=dega-1:-1:0
        temp=temp*eta+A(:,i*m+1:(i+1)*m);
    end    
    new_eta=max(eig(temp));
    if (new_eta < eta)
        eta_min=eta;
    else
        eta_max=eta;
    end
    eta=(eta_min+eta_max)/2;
end

if (nargout > 1)
    [V,D]=eig(temp');
    uT=V(:,find(sum(D,1)==max(sum(D,1))))';
end