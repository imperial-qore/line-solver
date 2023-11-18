function X=Q_Sylvest(U,T,B)
% X=Q_Sylvest(U,T,B) solves the equation X*kron(A,I)+BX=-I
% using a Hessenberg decomposition with kron(A,I)=U'*T*U

[V,NBAR]=hess(B);
%V'*B*V = NBAR

F=-V'*U;

n=size(F,2);

Y=[];
tempmat=zeros(n,n-1);
for k=1:n
    if (k==1)
        temp=F(:,k);
    else
        temp=F(:,k)-Y(:,1:k-1)*T(1:k-1,k);
    end
    Y(:,k)=(NBAR+eye(n)*T(k,k))\temp;
    % MATLAB checks that NBAR+T(k,k)*LBAR is an hessenberg
    % matrix and quickly reduces it to a triangular one which is
    % solved by backward substitution (TEST 6)
end
X=real(V*Y*U');