function G=MG1_EG(A,optVerbose)
%MG1_EG determines G directly if rank(A0)=1  

G=[];

m=size(A,1);
dega=size(A,2)/m-1;

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

if (drift < 1) % pos recurrent case
    if (rank(A(:,1:m))==1) % A0 = alpha * beta ?
        temp=find(sum(A(:,1:m),2)>0,1,'first'); % index of first nonzero row
        beta=A(temp,1:m)/sum(A(temp,1:m));
        G=ones(m,1)*beta;
    end
elseif (drift > 1) % transient case
    if (rank(A(:,1:m))==1) % A0 = alpha * beta ?
        if (optVerbose==1)
            Aold=A;
        end    
        for i=0:dega
            A(:,i*m+1:(i+1)*m)=diag(theta.^(-1))*A(:,i*m+1:(i+1)*m)'*diag(theta);
        end    
        etahat=GIM1_Caudal(A);
        temp=A(:,dega*m+1:end);
        for i=dega-1:-1:1
            temp=temp*etahat+A(:,i*m+1:(i+1)*m);
        end
        G=diag(theta.^(-1))*(A(:,1:m)*(eye(m)-temp)^(-1))'*...
            diag(theta);
        if (optVerbose==1)
            A=Aold;
        end    
    end
end

if (optVerbose==1)
    if (~isempty(G))
        Gcheck=A(:,dega*m+1:end);
        for j=dega-1:-1:0
            Gcheck=A(:,j*m+1:(j+1)*m)+Gcheck*G;
        end
        res_norm=norm(G-Gcheck,inf);
        fprintf('Final Residual Error for G: %d\n',res_norm);
    end
end
        