function G=MG1_RR(D,varargin)
%MG1_RR Ramaswami Reduction based Algorithm for M/G/1-Type Markov Chains 
%   [Bini,Meini,Ramaswami] 
%
%   G=MG1_RR(A) computes the minimal nonnegative solution to the 
%   matrix equation G = A0 + A1 G + A2 G^2 + A3 G^3 + ... + A_max G^max, 
%   where A = [A0 A1 A2 A3 ... A_max] has m rows and m*max columns and is
%   a nonnegative matrix, with (A0+A1+A2+...+A_max) irreducible and 
%   stochastic
%
%   Optional Parameters:
%   
%       MaxNumIt: Maximum number of iterations (default: 50)
%       Verbose: The residual error is printed at each step when set to 1,  
%                (default:0)
%       Mode: 'Direct' does not rely on the displacement structure, 
%             Requirements: memory O(m^2N^2), time O(m^3N^3)   
%             'DispStruct' makes use of the displacement structure,
%             Requirements: memory O(m^2N), time O(m^3N^2)
%             'DispStructFFT' uses the displacement structure and FFTs
%             Requirements: memory O(m^2N), time O(m^2NlogN+m^3N)
%             (default: 'DispStruct')

OptionNames=['Mode         '; 
             'MaxNumIt     ';
             'Verbose      '];
OptionTypes=['char   '; 
             'numeric';
             'numeric'];
OptionValues{1}=['Direct       '; 
                 'DispStruct   ';
                 'DispStructFFT'];
 
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end

% Default settings
%options.ProgressBar=0;
options.Mode='Direct';
options.MaxNumIt=50;
options.Verbose=0;

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

% check whether G is known explicitly
G=MG1_EG(D,options.Verbose);
if (~isempty(G))
    return
end

% Start RR
m=size(D,1);
N=size(D,2)/m-1;
D0=D(:,1:m);
vT=D(:,m+1:end);
uhat=zeros(N*m,m);
uhat(1:m,1:m)=eye(m);

if (strcmp(options.Mode,'Direct'))
    B=[zeros(m,N*m); eye((N-1)*m) zeros((N-1)*m,m)];
    check=min(norm(B,inf),norm(D0,inf));
    
    % ITERATIONS
    numit=0;
    while (check > 10^(-15) && numit < options.MaxNumIt)
        numit=numit+1;
        
        % step 1:
        S=inv(eye(m)-vT*uhat);
        ZZTuhat=[zeros(m,m) ; uhat(m+1:N*m,:)];
        T=vT*ZZTuhat;
        T=S*T;

        % step 2:
        newD0=D0*(S*vT(:,1:m)+T+eye(m))*D0;

        % step 3:
        % Use (7) of Theorem 3, Section 3.1
        newB=vT*B;
        newB=S*newB;
        newB=uhat*newB;
        newB=B+newB;
        newB=B*newB;
        % we continue with updating uhat and vT
        temp=uhat*(S*vT(:,1:m)+T+eye(m))*D0;
        newuhat=uhat+B*temp;
        temp=D0*S;
        temp=temp*vT;
        newvT=vT+temp*B;
        % update other variables
        B=newB;
        vT=newvT;
        uhat=newuhat;
        D0=newD0;
        
        check=min(norm(B,inf),norm(D0,inf));
        if (options.Verbose==1)
            fprintf('Check value of iteration %d equals %d \n',numit,check);
            drawnow;
        end
    end
else
    b=zeros(N*m,m);
    b(m+1:2*m,1:m)=eye(m);
    c1=zeros(m*N,m);
    c2=zeros(m*N,m);
    r1=zeros(m*N,m);
    r2=zeros(m*N,m);

    first=1;
    check=min(sum(sum(b)),norm(D0,inf));
    % ITERATIONS
    numit=0;
    while ((check > 10^(-15) || first == 1) && numit < options.MaxNumIt)
        numit=numit+1;
        first=0;

        % step 1:
        S=inv(eye(m)-vT*uhat);
        temp=[zeros(m,m) ; uhat(m+1:N*m,:)]; %% ZZTuhat
        T=vT*temp;
        T=S*T;

        % step 2:
        newD0=D0*(S*vT(:,1:m)+T+eye(m))*D0;

        % step 3:
        temp=uhat*(S*vT(:,1:m)+T+eye(m))*D0;
        if (strcmp(options.Mode,'DispStruct'))
            MG1_RR_Btemp;
        else
            MG1_RR_BtempFFT;
        end
        newuhat=uhat+temp;
        temp=D0*S;
        temp=temp*vT;
        if (strcmp(options.Mode,'DispStruct'))
            MG1_RR_tempB;
        else
            MG1_RR_tempBFFT;
        end
        newvT=vT+temp;

        % step 4:
        % we start by calculating (I-A)^(-1) b using the structure of (I-A)^-1
        temp=[vT*b; b(1:m,1:m)];
        temp=[S T; S eye(m)+T]*temp;
        newb(1:m,1:m)=b(1:m,1:m)+temp(1:m,:);
        newb(m+1:N*m,1:m)=b(m+1:end,1:m)+uhat(m+1:end,:)*temp(m+1:end,:);
        temp=newb;
        if (strcmp(options.Mode,'DispStruct'))
            MG1_RR_Btemp;
        else
            MG1_RR_BtempFFT;
        end
        newb=temp;

        % step 5, part 1:
        % H
        e1pZu=[eye(m); uhat(m+1:end,:)];
        Z2u=[zeros(2*m,m); uhat(m+1:(N-1)*m,:)];
        e2pZ2u=[zeros(m,m); eye(m); uhat(m+1:(N-1)*m,:)];
        H=[e1pZu e2pZ2u*S e2pZ2u*T+Z2u];

        clear e1pZu;
        clear Z2u;
        clear e2pZ2u;

        % KT
        temp=[vT(:,m+1:end) zeros(m,m)];
        KT=[S*temp; -vT; [-eye(m) zeros(m,(N-1)*m)]];

        % W
        W=[c1 c2];
        if (strcmp(options.Mode,'DispStruct'))
            temp=H;
            MG1_RR_Btemp;
            W(1:N*m,2*m+1:5*m)=temp;
        else
            for colp=0:2
                temp=H(:,colp*m+1:(colp+1)*m);
                MG1_RR_BtempFFT;
                W(1:N*m,(2+colp)*m+1:(3+colp)*m)=temp;
            end
        end
        clear H;
        temp=[vT*[c1 c2]; [c1(1:m,1:m) c2(1:m,1:m)]];
        temp=[S T; S eye(m)+T]*temp;
        temp2(1:m,1:2*m)=[c1(1:m,1:m) c2(1:m,1:m)]+temp(1:m,:);
        temp2(m+1:N*m,1:2*m)=[c1(m+1:end,1:m) c2(m+1:end,1:m)]+...
            uhat(m+1:end,:)*temp(m+1:end,:);
        if (strcmp(options.Mode,'DispStruct'))
            temp=temp2;
            clear temp2;
            MG1_RR_Btemp;
            W(1:N*m,5*m+1:7*m)=temp;
        else
            supertemp2=temp2;
            clear temp2;
            for colp=0:1
                temp=supertemp2(:,colp*m+1:(1+colp)*m);
                MG1_RR_BtempFFT;
                W(1:N*m,(5+colp)*m+1:(6+colp)*m)=temp;
            end
            clear supertemp2;
        end
        clear temp;

        % Step 6, part 1 (this order allows us to clear W more rapidly)
        [Q1,R1]=qr(W,0);
        clear W;

        % YT
        YT(5*m+1:7*m,1:N*m)=[r1 r2]';
        if (strcmp(options.Mode,'DispStruct'))
            temp=KT;
            MG1_RR_tempB;
            YT(2*m+1:5*m,1:N*m)=temp;
        else
            for rowp=0:2
                temp=KT(rowp*m+1:(1+rowp)*m,:);
                MG1_RR_tempBFFT;
                YT((2+rowp)*m+1:(3+rowp)*m,1:N*m)=temp;
            end
        end    
        clear KT;
        Zu=[zeros(m,m); uhat(m+1:end,:)];
        temp=[[r1(1:m,1:m) r2(1:m,1:m)]' [r1 r2]'*Zu];
        temp=temp*[S T; S eye(m)+T];
        temp2=temp(:,1:m)*vT;
        temp2(:,1:m)=temp2(:,1:m)+temp(:,m+1:2*m);
        if (strcmp(options.Mode,'DispStruct'))
            temp=temp2+[r1 r2]';
            clear temp2;
            MG1_RR_tempB;
            YT(1:2*m,1:N*m)=temp;
        else
            supertemp2=temp2+[r1 r2]';
            clear temp2;
            for rowp=0:1
                temp=supertemp2(rowp*m+1:(1+rowp)*m,:);
                MG1_RR_tempBFFT;
                YT(rowp*m+1:(1+rowp)*m,1:N*m)=temp;
            end
            clear supertemp2;
        end    
        clear temp;

        % Step 6, part 2:
        [Q2,R2]=qr(YT',0);
        clear YT;

        % Step 7:
        [U,Xsi,V]=svd(R1*R2');
        for temp=1:2*m
            first2m(temp)=Xsi(temp,temp);
        end
        clear Xsi;
        clear R1;
        clear R2;

        % Step 8:
        temp=Q1*U(:,1:2*m)*diag(first2m);
        clear first2m;
        clear Q1;
        clear U;
        c1=temp(:,1:m);
        c2=temp(:,m+1:2*m);
        temp=(V(:,1:2*m)'*Q2')';
        clear Q2;
        clear V;
        r1=temp(:,1:m);
        r2=temp(:,m+1:2*m);

        % update other variables
        b=newb;
        vT=newvT;
        uhat=newuhat;
        D0=newD0;

        clear newb;
        clear newvT;
        clear newuhat;
        clear newD0;

        check=min(sum(sum(b)),norm(D0,inf));
        if (options.Verbose==1)
            fprintf('Check value of iteration %d equals %d \n',numit,check);
            drawnow;
        end
        
    end
end
if (numit == options.MaxNumIt && check > 10^(-15))
    warning('Maximum Number of Iterations %d reached',numit);
end

% Compute G by means of formula at the end of Section 3.
temp=eye(m)-D(:,m+1:end)*uhat;
G=temp^(-1)*D(:,1:m);

if (options.Verbose==1)
    temp=D(:,end-m+1:end);
    for i=size(D,2)/m-1:-1:1
        temp=D(:,(i-1)*m+1:i*m)+temp*G;
    end
    res_norm=norm(G-temp,inf);
    fprintf('Final Residual Error for G: %d\n',res_norm);
end

 

    
