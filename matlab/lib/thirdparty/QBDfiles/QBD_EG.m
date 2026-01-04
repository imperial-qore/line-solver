function [G,R,U]=QBD_EG(A0,A1,A2,optVerbose,nargout_caller)
%QBD_EG determines G directly if rank(A0)=1 or rank(A2)=1 

G=[];
R=[];
U=[];

m=size(A1,1);
theta=stat(A0+A1+A2);
drift=theta*sum(A0,2)-theta*sum(A2,2);
if (drift > 0) % pos recurrent case
    if (rank(A0)==1) % A0 = alpha * beta ?
        temp=min(find(sum(A0,2)>0)); % find index of first nonzero row
        beta=A0(temp,:)/sum(A0(temp,:));
        G=ones(m,1)*beta;
        if (nargout_caller > 1)
            R=A2*(eye(m)-(A1+A2*G))^(-1);
        end    
    elseif (rank(A2)==1) % A2 = alpha * beta ?
        eta=QBD_Caudal(A0,A1,A2);
        R=A2*(eye(m)-A1-eta*A0)^(-1);
        G=(eye(m)-(A1+R*A0))^(-1)*A0;   
    end    
elseif (drift < 0) % transient case
    if (rank(A2)==1) % A2 = alpha * beta ?
        alpha=A2*ones(m,1);
        R=alpha*theta/(theta*alpha);
        G=(eye(m)-(A1+R*A0))^(-1)*A0;   
    elseif (rank(A0)==1) % A0 = alpha * beta ?
        A0hat=diag(theta.^(-1))*A2'*diag(theta);
        A1hat=diag(theta.^(-1))*A1'*diag(theta);
        A2hat=diag(theta.^(-1))*A0'*diag(theta);
        etahat=QBD_Caudal(A0hat,A1hat,A2hat);
        G=diag(theta.^(-1))*(A2hat*(eye(m)-A1hat-etahat*A0hat)^(-1))'*...
            diag(theta);
        if (nargout_caller > 1)
            R=A2*(eye(m)-(A1+A2*G))^(-1);
        end
    end
end

if (~isempty(R) & nargout_caller > 2)
    U=A1+R*A0;
end

if (optVerbose==1)
    if (~isempty(G))
        res_norm=norm(G-A0-(A1+A2*G)*G,inf);
        fprintf('Final Residual Error for G: %d\n',res_norm);
    end
    if (~isempty(R) & nargout_caller > 1)
        res_norm=norm(R-A2-R*(A1+R*A0),inf); 
        fprintf('Final Residual Error for R: %d\n',res_norm);
    end
    if (~isempty(U))
        res_norm=norm(U-A1-A2*(eye(m)-U)^(-1)*A0,inf);
        fprintf('Final Residual Error for U: %d\n',res_norm);
    end
end    
        