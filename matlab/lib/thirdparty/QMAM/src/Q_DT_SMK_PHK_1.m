function [ql,qlT,wait,waitT,soj,sojT]=Q_DT_SMK_PHK_1(D,alpha,S,varargin)
%[ql,qlT,wait,waitT,soj,sojT]=Q_DT_SMK_PHK_1(D,alpha,S) computes the 
%   Queuelength, Waiting Time and Sojourn time distribution of an
%   SM[K]/PH[K]/1/FCFS queue (per type and overall)
%   
%   INPUT PARAMETERS:
%   * SM[K] arrival process (with m states)
%     D{k}=[Dk_0 Dk_1 ... Dk_Imax], for k=1...K, where the entry (j,j') of
%     the mxm matrix Dk_i, holds the probabilities of having a type k
%     arrival, with an interarrival time of i time slots, while the 
%     underlying phase changes from j to j'
%
%   * PH[K] services time distributions (i=1:K)
%     alpha{k} holds the 1xmk alpha vector of the PH service time of
%     type k customers S{k} holds the mixmi matrix S of the 
%     PH service time of a type k customer
%
%   RETURN VALUES:
%   * Per type queue length distribution, 
%     ql{k}(i) = Prob[(i-1) type k customers in the queue]
%   * Overall queue length distribution,
%     qlT(i) = Prob[(i-1) customers in the queue (any type)]
%   * Per type waiting time distribution, 
%     wait{k}(i) = Prob[a type k customers has waiting time = (i-1)]
%   * Overall waiting time distribution,
%     waitT(i) = Prob[a customer (of any type) has wainting time = (i-1)]
%   * Per type sojourn time distribution, 
%     soj{k}(i) = Prob[a type k customers has sojourn time = (i-1)]
%   * Overall sojourn time distribution,
%     sojT(i) = Prob[a customer (of any type) has sojourn time = (i-1)]
%
%
%   OPTIONAL PARAMETERS:
%
%       Mode: The underlying function to compute the R matrix of the 
%             underlying GIM1-type Markov chain can be selected using 
%             the following parameter values (default: 'CR')
%               
%               'CR' : Cyclic Reduction [Bini, Meini]
%               'FI' : Functional Iterations [Neuts]
%               'IS' : Invariant Subspace [Akar, Sohraby]
%               'RR' : Ramaswami Reduction [Bini, Meini, Ramaswami]
%
%       MaxNumComp: Maximum number of components for the vectors containig
%           the performance measure.
%       
%       Verbose: When set to 1, the progress of the computation is printed
%           (default:0).
%
%       Optfname: Optional parameters for the underlying function fname.
%           These parameters are included in a cell with one entry holding
%           the name of the parameter and the next entry the parameter
%           value. In this function, fname can be equal to:
%               'GIM1_R' : Options for the GIM1_R functions
%
%   USES: GIM1_R of the SMCSolver tool

OptionNames=[
             'Mode              ';
             'MaxNumComp        ';
             'Verbose           ';
             'OptGIM1_R         '];
OptionTypes=[
             'char   ';
             'numeric';
             'numeric';
             'cell   '];
         
OptionValues{1}=['CR';
                 'FI';
                 'IS';
                 'RR'];
             
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.Mode='CR';
options.MaxNumComp = 1000;
options.Verbose = 0;
options.OptGIM1_R=cell(0);


% Parse Parameters
Q_DT_SMK_ParsePara(D,'D',1);
K=size(alpha,2);
for i=1:K
    Q_DT_PH_ParsePara(alpha{i}, ['alpha_' int2str(i)], S{i}, ['S_', int2str(i)] )
end

% Parse Optional Parameters
options=Q_ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

% Determine constants
m=size(D{1},1);
K=size(alpha,2);
smk=zeros(1,K+1);
smk(1)=0;
for i=1:K
    smk(i+1)=smk(i)+size(alpha{i},2);
end    
mser=smk(K+1);
mtot=mser*m;

% Test the load of the queue
Dsum=zeros(m,m);
meanI=zeros(m,m);
maxI=1;
maxIs=ones(1,K);
for k=1:K
    Dsums{k}=zeros(m,m);
    for i=1:size(D{k},2)/m 
        Dsums{k}=Dsums{k}+D{k}(:,(i-1)*m+1:i*m);
        meanI=meanI+i*D{k}(:,(i-1)*m+1:i*m);
    end
    Dsum=Dsum+Dsums{k};
    maxIs(k)=size(D{k},2)/m;
end    
maxI=max(maxIs);
pi=stat(Dsum);
lambda=1/(pi*meanI*ones(m,1));
lambdas=zeros(1,K);
mus=zeros(1,K);
for k=1:K
   lambdas(k)=lambda*pi*Dsums{k}*ones(m,1);
   means(k)=sum(alpha{k}*(eye(size(alpha{k},2))-S{k})^(-1));
end   
load=lambdas*means';
if load >= 1
    error('MATLAB:Q_DT_SMK_PHK_1:LoadExceedsOne',...
                        'The load %d of the system exceeds one',load);
end    

% Construct building blocks A0, A1, ... AImax
Dmats=zeros(m*maxI,m*K);
Dmatssum=zeros(m*maxI,m);
for k=1:k
    for i=1:size(D{k},2)/m 
        Dmats((i-1)*m+1:i*m,(k-1)*m+1:k*m)=D{k}(:,(i-1)*m+1:i*m);
        Dmatssum((i-1)*m+1:i*m,:)=Dmatssum((i-1)*m+1:i*m,:)+D{k}(:,(i-1)*m+1:i*m);
    end
end

Tser=zeros(mser,mser);
tser=zeros(mser,1);
for k=1:K
    Tser(smk(k)+1:smk(k+1),smk(k)+1:smk(k+1))=S{k};
    tser(smk(k)+1:smk(k+1),1)=1-sum(S{k},2);
end

A0=kron(Tser,eye(m));
Temp=zeros(m*maxI,mtot); % i > 0
for k=1:K
    Temp(:,m*smk(k)+1:m*smk(k+1))=kron(alpha{k},Dmats(:,(k-1)*m+1:k*m));
end
Ais=zeros(mtot*maxI,mtot); % i > 0
for i=1:maxI
    Ais((i-1)*mtot+1:i*mtot,:)=kron(tser,Temp((i-1)*m+1:i*m,:));
end

%Prepare input GIM1_R
A=zeros(mtot,mtot*(maxI+1));
A(:,1:mtot)=A0;
for i=1:maxI
    A(:,i*mtot+1:(i+1)*mtot)=Ais((i-1)*mtot+1:i*mtot,:);
end    

switch options.Mode
    case 'CR'
        R=GIM1_R(A,'CR',options.OptGIM1_R{:});
    case 'FI'
        R=GIM1_R(A,'FI',options.OptGIM1_R{:});
    case 'IS'
        R=GIM1_R(A,'IS',options.OptGIM1_R{:});
    case 'RR'
        R=GIM1_R(A,'RR',options.OptGIM1_R{:});
end

% compute theta_tot
Asum=A0;
for i=1:maxI
    Asum=Asum+A(:,i*mtot+1:(i+1)*mtot);
end    
theta_tot=stat(Asum);

% compute pi's with preallocation based on Caudal characteristic
pi=zeros(ceil(log(10^(-12))/log(max(eig(R)))),mtot);
pi(1,:)=load*theta_tot*(eye(mtot)-R);
sumpi=sum(pi(1,:));
i=2;
while (sumpi < load*(1-10^(-12)) && i < 1+options.MaxNumComp )
    pi(i,1:mtot)=pi(i-1,1:mtot)*R;
    sumpi=sumpi+sum(pi(i,1:mtot));
    i=i+1;
    if (~mod(i,options.Verbose))
        fprintf('Accumulated mass after %d iterations: %d\n',i,sumpi);
        drawnow;
    end
end    
if (i == 1+options.MaxNumComp)
    warning('Maximum Number of Components %d reached',i-1);
end

% compute sojourn times per type
piS=sum(reshape(pi',m,size(pi,1)*mser),1);
piS=reshape(piS,mser,size(piS,2)/mser);
for k=1:K
    soj{k}=zeros(1,size(piS,2)+1);
end
for i=1:size(piS,2)
    for k=1:K %order: good for caching (MATLAB is columnwise)
        soj{k}(i+1)=piS(smk(k)+1:smk(k+1),i)'*tser(smk(k)+1:smk(k+1),1);
    end
end 
for k=1:K
    soj{k}=soj{k}/lambdas(k);
end

% compute overall sojourn time 
sojT=zeros(size(soj{1}));
for i=1:K
   sojT=sojT+lambdas(i)*soj{i}; 
end
sojT=sojT/lambda;

% compute the waiting time per type (no deconvolution -> not stable)
maxwait=size(pi,1)-1;
piT=pi*kron(tser,eye(m));
piT=reshape(piT',1,m*size(pi,1));
Dmatse=zeros(maxI*m,K);
for k=1:K
    Dmatse(:,k)=sum(Dmats(:,(k-1)*m+1:k*m),2);
end
for k=1:K
    wait{k}=zeros(1,maxwait+1);
    for n=1:maxwait
        wait{k}(n+1)=piT(n*m+1:min(maxwait+1,n+maxIs(k))*m)*...
        Dmatse(1:m*min(maxIs(k),maxwait+1-n),k);
    end
    wait{k}=wait{k}/lambdas(k);
    wait{k}(1)=1-sum(wait{k});
end

%compute overall waiting time
waitT=zeros(size(wait{1}));
for i=1:K
   waitT(1,1:size(wait{i},2))=waitT(1,1:size(wait{i},2))+lambdas(i)*wait{i}; 
end
waitT=waitT/sum(lambdas);

% compute queue length per type (new arrivals at current time are not included) 
piTA=zeros(size(pi,1),m*K);
for k=1:K
    for j=1:m
        piTA(:,m*(k-1)+j)=sum(pi(:,m*smk(k)+j:m:m*smk(k+1)),2);
        % elements piTA(i,m*(k-1)+1:m*k) give type k in service
    end    
end
piA=zeros(size(pi,1),m);
for j=1:m
    piA(:,j)=sum(piTA(:,j:m:end),2); 
end 
for k=1:K
    Dtemps=zeros(m*maxI,size(pi,1));
    piAk=piA-piTA(:,m*(k-1)+1:m*k); % no type k in service
    ql{k}=zeros(1,size(pi,1)+1);
    ql{k}(1)=1-load;
    
    Dkmats=zeros(m,m*maxIs(k));
    for i=1:maxIs(k)
        Dkmats(:,(maxIs(k)-i)*m+1:(maxIs(k)-i+1)*m)=...
            Dmats((i-1)*m+1:i*m,(k-1)*m+1:k*m);
    end    % Dkmats = [Dk(maxIs(k)) ... Dk(2) Dk(1)]
    Dnokmats=zeros(m,m*maxI);
    for i=1:maxI
        Dnokmats(:,(maxI-i)*m+1:(maxI-i+1)*m)=...
            Dmatssum((i-1)*m+1:i*m,:)-Dmats((i-1)*m+1:i*m,(k-1)*m+1:k*m);
    end    % Dnokmats = [Dnok(maxI) ... Dnok(2) Dnok(1)]
    
    Dtemps(1:m,1)=ones(m,1);
    for i=1:maxIs(k)-1
        ql{k}(1:i)=ql{k}(1:i)+piAk(i,:)*Dtemps((i-1)*m+1:i*m,1:i);
        % type k in service
        ql{k}(2:i+1)=ql{k}(2:i+1)+piTA(i,m*(k-1)+1:m*k)*...
            Dtemps((i-1)*m+1:i*m,1:i);
        % update Dtemps
        Dtemps(i*m+1:(i+1)*m,1:(i+1))=...
            [zeros(m,1) Dkmats(:,m*(maxIs(k)-i)+1:m*maxIs(k))*Dtemps(1:m*i,1:i)]+...
            [Dnokmats(:,m*(maxI-i)+1:m*maxI)*Dtemps(1:m*i,1:i) zeros(m,1)];
        Dtemps(i*m+1:(i+1)*m,1)=ones(m,1)-sum(Dtemps(i*m+1:(i+1)*m,2:(i+1)),2);
    end    
    for i=maxIs(k):maxI-1
        ql{k}(1:i)=ql{k}(1:i)+piAk(i,:)*Dtemps((i-1)*m+1:i*m,1:i);
        % type k in service
        ql{k}(2:i+1)=ql{k}(2:i+1)+piTA(i,m*(k-1)+1:m*k)*...
            Dtemps((i-1)*m+1:i*m,1:i);
        % update Dtemps
        Dtemps(i*m+1:(i+1)*m,1:(i+1))=...
            [zeros(m,1) Dkmats*Dtemps((i-maxIs(k))*m+1:m*i,1:i)]+...
            [Dnokmats(:,m*(maxI-i)+1:m*maxI)*Dtemps(1:m*i,1:i) zeros(m,1)];
        Dtemps(i*m+1:(i+1)*m,1)=ones(m,1)-sum(Dtemps(i*m+1:(i+1)*m,2:(i+1)),2);
    end
    for i=maxI:size(pi,1)
        ql{k}(1:i)=ql{k}(1:i)+piAk(i,:)*Dtemps((maxI-1)*m+1:maxI*m,1:i);
        % type k in service
        ql{k}(2:i+1)=ql{k}(2:i+1)+piTA(i,m*(k-1)+1:m*k)*...
            Dtemps((maxI-1)*m+1:maxI*m,1:i);
        % update Dtemps
        Dtemps(maxI*m+1:(maxI+1)*m,1:(i+1))=...
            [zeros(m,1) Dkmats*Dtemps((maxI-maxIs(k))*m+1:m*maxI,1:i)]+...
            [Dnokmats*Dtemps(:,1:i) zeros(m,1)];
        Dtemps(maxI*m+1:(maxI+1)*m,1)=...
            ones(m,1)-sum(Dtemps(maxI*m+1:(maxI+1)*m,2:(i+1)),2);
        Dtemps(1:m,:)=[];
    end
    ql{k}=ql{k}(1:max(find(ql{k}>10^(-14))));
end

% compute overall queue length
Dtemps=zeros(m*maxI,size(pi,1));
qlT=zeros(1,size(pi,1)+1);
qlT(1)=1-load;

DTmats=zeros(m,m*maxI);
for i=1:maxI
    DTmats(:,(maxI-i)*m+1:(maxI-i+1)*m)=Dmatssum((i-1)*m+1:i*m,:);
end

Dtemps(1:m,1)=ones(m,1);
for i=1:maxI-1
    % customer in service
    qlT(2:i+1)=qlT(2:i+1)+piA(i,:)*Dtemps((i-1)*m+1:i*m,1:i);
    % update Dtemps
    Dtemps(i*m+1:(i+1)*m,1:(i+1))=...
        [zeros(m,1) DTmats(:,m*(maxI-i)+1:m*maxI)*Dtemps(1:m*i,1:i)];
    Dtemps(i*m+1:(i+1)*m,1)=ones(m,1)-sum(Dtemps(i*m+1:(i+1)*m,2:(i+1)),2);
end
for i=maxI:size(pi,1)
    % customer in service
    qlT(2:i+1)=qlT(2:i+1)+piA(i,:)*Dtemps((maxI-1)*m+1:maxI*m,1:i);
    % update Dtemps
    Dtemps(maxI*m+1:(maxI+1)*m,1:(i+1))=[zeros(m,1) DTmats*Dtemps(:,1:i)];
    Dtemps(maxI*m+1:(maxI+1)*m,1)=...
        ones(m,1)-sum(Dtemps(maxI*m+1:(maxI+1)*m,2:(i+1)),2);
    Dtemps(1:m,:)=[];
end
qlT=qlT(1:max(find(qlT>10^(-14))));
