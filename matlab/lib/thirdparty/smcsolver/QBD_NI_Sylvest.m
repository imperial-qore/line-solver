function X=QBD_NI_Sylvest(A,B,C,D)
% this function solves the equation AXB+CX=D
% using a complex schur and a hessenberg-triangular decomposition

 [LBAR,NBAR,W,V]=hess(A,C);
 %W*A*V = LBAR, W*C*V = NBAR
 [U,T]=schur(B,'complex');
 %U'*T*U = B
 F=W*D*U;
 
 n=size(F,2);
 
 Y=[];
 tempmat=zeros(n,n-1);
 for k=1:n
     if (k==1)
         temp=F(:,k);
     else    
         tempmat(:,k-1)=-LBAR*Y(:,k-1);
         temp=F(:,k)+sum(tempmat(:,1:k-1).*kron(ones(n,1),T(1:k-1,k).'),2);
     end
     Y(:,k)=(NBAR+T(k,k)*LBAR)\temp;
     % MATLAB checks that NBAR+T(k,k)*LBAR is an hessenberg
     % matrix and quickly reduces it to a triangular one which is 
     % solved by backward substitution (TEST 6)
 end
 X=real(V*Y*U');