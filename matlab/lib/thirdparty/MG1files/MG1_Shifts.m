function [hatA,drift,tau,v]=MG1_Shifts(A,ShiftType)
% this function performs the shift operation for the MG1 case
% input: A = [A0 A1 A2 ... A_maxd]
% includes: one, tau and dbl

m=size(A,1);
v=zeros(m,1); % default value if redundant
tau=1; % default value if redundant
maxd=size(A,2)/m-1;
sumA=A(:,maxd*m+1:end);
beta=sum(sumA,2);
% beta = (A_maxd)e + (A_maxd + A_maxd-1)e + ... + (Amaxd+...+A1)e
for i=maxd-1:-1:1
    sumA=sumA+A(:,i*m+1:(i+1)*m);
    beta=beta+sum(sumA,2);
end
sumA=sumA+A(:,1:m);
theta=stat(sumA);
drift=theta*beta;

% if drift < 1 : positive recurrent
hatA=zeros(m,m*(maxd+1));
if (drift < 1)
    if (strcmp(ShiftType,'tau')|strcmp(ShiftType,'dbl'))  % shift tau to infinity
            [tau,uT]=MG1_Decay(A);
            A(:,m+1:2*m)=A(:,m+1:2*m)-eye(m);
            uT=uT/sum(uT);
            rowhatA=zeros(1,m*(maxd+1));
            rowhatA(1,maxd*i:end)=uT*A(:,maxd*i:end);
            for i=maxd-1:-1:0
                rowhatA(1,i*m+1:(i+1)*m)=tau*rowhatA(1,(i+1)*m+1:(i+2)*m)+...
                    uT*A(:,i*m+1:(i+1)*m); 
            end
            hatA=A-ones(m,1)*rowhatA;
    end       
    if (strcmp(ShiftType,'dbl')) % shift one to zero
        A=hatA;
        % e is also the right eigenvector of hatA(1)
        % as the shift-tau does not influence G and Ge=e,
        % implying that hatA(1)e = e as G=hatA(G)
    end
    if (strcmp(ShiftType,'one')) % shift one to zero
            A(:,m+1:2*m)=A(:,m+1:2*m)-eye(m);
    end
    if (strcmp(ShiftType,'one')|strcmp(ShiftType,'dbl')) % shift one ot zero
            colhatA=zeros(m,maxd+1);
            colhatA(:,1)=sum(A(:,1:m),2); % colhatA(:,1) = (A0)e
            for i=1:maxd
                colhatA(:,i+1)=colhatA(:,i)+sum(A(:,i*m+1:(i+1)*m),2);
                % colhatA(:,i+1) = (A0+A1+...+Ai)e
            end
            hatA=A-kron(colhatA,ones(1,m)/m); % hatAi = Ai - (A0+A1+...+Ai)e*uT
    end
else
    if (strcmp(ShiftType,'one')|strcmp(ShiftType,'dbl')) % shift one to infinity
        A(:,m+1:2*m)=A(:,m+1:2*m)-eye(m);
        rowhatA=zeros(1,m*(maxd+1));
        rowhatA(1,maxd*i:end)=theta*A(:,maxd*i:end);
        for i=maxd-1:-1:0
            rowhatA(1,i*m+1:(i+1)*m)=rowhatA(1,(i+1)*m+1:(i+2)*m)+...
                theta*A(:,i*m+1:(i+1)*m); % rowhatAi = theta(Amaxd+...+Ai)
        end
        hatA=A-ones(m,1)*rowhatA;
    end  
    if (strcmp(ShiftType,'dbl')) % shift one to infinity
        A=hatA;
        A(:,m+1:2*m)=A(:,m+1:2*m)+eye(m);
        % v is also the right eigenvector of hatA(tau)
        % as the shift-one does not influence G and Gv=tau*v,
        % implying that hatA(tau)v = tau*v as G=hatA(G)
    end
    if (strcmp(ShiftType,'tau')|strcmp(ShiftType,'dbl'))  % shift tau to zero
        [tau,v]=GIM1_Caudal(A);
        A(:,m+1:2*m)=A(:,m+1:2*m)-eye(m);
        v=v/sum(v);
        colhatA=zeros(m,maxd+1);
        colhatA(:,1)=A(:,1:m)*v; % colhatA(:,1) = (A0)v
        for i=1:maxd
            colhatA(:,i+1)=colhatA(:,i)*tau^(-1)+A(:,i*m+1:(i+1)*m)*v;
        end
        hatA=A-kron(colhatA,ones(1,m)); 
    end      
end    
hatA(:,m+1:2*m)=hatA(:,m+1:2*m)+eye(m);

