function eta=GIM1_CAUDAL(A)
%GIM1_CAUDAL Computes the Spectral Radius of R 
%
%   eta=GIM1_CAUDAL(A) computes the dominant eigenvalue of the
%   matrix R, the smallest nonnegative solution to 
%   R= A0 + R A1 + R^2 A2 + ... + R^max Amax

m=size(A,1);
dega=size(A,2)/m-1;

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
