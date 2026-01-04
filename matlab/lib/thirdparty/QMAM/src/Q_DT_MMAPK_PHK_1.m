function [ql,qlT,wait,waitT,soj,sojT]=Q_DT_MMAPK_PHK_1(D0,D,alpha,S,varargin)
%[ql,qlT,wait,waitT,soj,sojT]=Q_DT_MMAPK_PHK_1(D0,D,alpha,S) computes the 
%   Queuelength, Waiting time and Sojourn time distribution of a 
%   Discrete-Time MMAP[K]/PH[K]/1/FCFS queue (per type and overall)
%   
%   INPUT PARAMETERS:
%   * MMAP[K] arrival process (with m states)
%     D{i} holds the mxm matrix D_i, together with the mxm matrix D0 
%     these characterize the MMAP[K] arrival process
%
%   * PH[K] services time distributions (i=1:K)
%     alpha{i} holds the 1xmi alpha vector of the PH service time of
%     type i customers
%     S{i} holds the mixmi matrix S of the PH service time of a type
%     i customer
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
%   OPTIONAL PARAMETERS:
%       Mode: The underlying function to compute the R matrix of the 
%             underlying QBD can be selected using the following 
%             parameter values (default: 'CR')
%               'CR' : Cyclic Reduction [Bini, Meini]
%               'FI' : Functional Iterations [Neuts]
%               'IS' : Invariant Subspace [Akar, Sohraby]
%               'LR' : Logaritmic Reduction [Latouche, Ramaswami]
%               'NI' : Newton Iteration
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
%               'QBD_CR' : Options for Cyclic Reduction [Bini, Meini]
%               'QBD_FI' : Options for Functional Iterations [Neuts]
%               'QBD_IS' : Options for Invariant Subspace [Akar, Sohraby]
%               'QBD_LR' : Options for Logaritmic Reduction [Latouche, Ramaswami]
%               'QBD_NI' : Options for Newton Iteration
%
%   USES: QBD Solver and QBD_pi of the SMCSolver tool


OptionNames=[
             'Mode              ';   
             'MaxNumComp        ';
             'Verbose           ';
             'OptQBD_CR         '; 
             'OptQBD_FI         '; 
             'OptQBD_IS         '; 
             'OptQBD_LR         '; 
             'OptQBD_NI         '];
OptionTypes=[
             'char   ';
             'numeric';
             'numeric';
             'cell   '; 
             'cell   '; 
             'cell   '; 
             'cell   '; 
             'cell   '];
         
OptionValues{1}=['CR';
                 'FI';
                 'IS';
                 'LR';
                 'NI'];
             
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.Mode='CR';
options.MaxNumComp = 1000;
options.Verbose = 0;
options.OptQBD_CR=cell(0);
options.OptQBD_FI=cell(0);
options.OptQBD_IS=cell(0);
options.OptQBD_LR=cell(0);
options.OptQBD_NI=cell(0);

% Parse Parameters
K=size(alpha,2);
Q_DT_MMAPK_ParsePara(D0,'D0',D,'D');
for i=1:K
    Q_DT_PH_ParsePara(alpha{i}, ['alpha_' int2str(i)], S{i}, ['S_', int2str(i)] )
end

% Parse Optional Parameters
options=Q_ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

% Check Unused Optional Parameters
Q_CheckUnusedParaQBD(options);

% Determine constants
m=size(D0,1);
K=size(alpha,2);
size(alpha{1},2);
smk=zeros(1,K+1);
smk(1)=0;
for i=1:K
    smk(i+1)=smk(i)+size(alpha{i},2);
end    
mser=smk(K+1);
mtot=mser*m;

% Test the load of the queue
Dsum=D0;
for i=1:K
    Dsum=Dsum+D{i};
end    
pi=stat(Dsum);
lambdas=zeros(1,K);
mus=zeros(1,K);
for i=1:K
   lambdas(i)=pi*D{i}*ones(m,1);
   means(i)=sum(alpha{i}*(eye(size(alpha{i},2))-S{i})^(-1));
end   
lambdas
means
lambdas.*means
load=lambdas*means'
if load >= 1
    error('MATLAB:Q_DT_MMAPK_PHK_1:LoadExceedsOne',...
                        'The load %d of the system exceeds one',load);
end    

% Construct building blocks
L=zeros(m,mtot);
Tser=zeros(mser,mser);
tser=zeros(mser,1);

for i=1:K
    L(:,m*smk(i)+1:m*smk(i+1))=kron(alpha{i},D{i});
    Tser(smk(i)+1:smk(i+1),smk(i)+1:smk(i+1))=S{i};
    tser(smk(i)+1:smk(i+1),1)=1-sum(S{i},2);
end

% Compute QBD blocks A0, A1 and A2
A0=[zeros(m,m+mtot); zeros(mtot,m) kron(Tser,eye(m))];
A1=[zeros(m,m) L; zeros(mtot,m) kron(tser,L)];
A2=[D0 zeros(m,mtot); kron(tser,D0) zeros(mtot,mtot)];

B0=[zeros(m,m) L];
B1=D0;
B2=[D0; kron(tser,D0)];

switch options.Mode
    case 'CR'
        [G,R]=QBD_CR(A2,A1,A0,options.OptQBD_CR{:});
    case 'FI'
        [G,R]=QBD_FI(A2,A1,A0,options.OptQBD_FI{:});
    case 'LR'
        [G,R]=QBD_LR(A2,A1,A0,options.OptQBD_LR{:});
    case 'IS'
        [G,R]=QBD_IS(A2,A1,A0,options.OptQBD_IS{:});
    case 'NI'
        [G,R]=QBD_NI(A2,A1,A0,options.OptQBD_NI{:});
end
pi=QBD_pi(B2,B1,R,'MaxNumComp',options.MaxNumComp,'Verbose', options.Verbose,'Boundary',[B0; A1+R*A2]);

% remove additional states
pi0=pi(1:m);
pi=reshape(pi(m+1:end),mtot+m,size(pi(m+1:end),2)/(mtot+m));
c=sum(sum(pi(1:m,:)));
pi0=pi0/(1-c);
pi=pi(m+1:end,:)/(1-c);

% compute sojourn times per type
piS=sum(reshape(pi,m,size(pi,2)*mser),1);
piS=reshape(piS,mser,size(piS,2)/mser);
for k=1:K
    soj{k}=zeros(1,size(piS,2)+1);
end
for i=1:size(pi,2)
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
sojT=sojT/sum(lambdas);

% compute the waiting time per type
piW=pi'*kron(tser,eye(m));
maxwait=size(piW,1)-1;
piW=reshape(piW',1,m*size(pi,2));
Dmatse=zeros((maxwait+1)*m,K);
for k=1:K
    temp=sum(D{k},2);
    for i=1:maxwait+1
        Dmatse((i-1)*m+1:i*m,k)=temp;
        temp=D0*temp;
    end    
end
for k=1:K
    wait{k}=zeros(1,maxwait+1);
    for n=1:maxwait
        wait{k}(n+1)=piW(n*m+1:(maxwait+1)*m)*...
        Dmatse(1:m*(maxwait+1-n),k);
    end
    wait{k}=wait{k}/lambdas(k);
    wait{k}(1)=1-sum(wait{k});
end

% old waiting time -> sometimes unstable
% for i=1:K
%     temp=alpha{i};
%     prob=0;
%     sprob=prob;
%     j=1;
%     while (sprob < 1-10^(-10))
%         j=j+1;
%         prob(j)=temp*tser(smk(i)+1:smk(i+1));
%         temp=temp*S{i};
%         sprob=sprob+prob(j);
%     end
%     prob=prob(min(find(prob>10^(-12))):end);
%     stemp=soj{i}(min(find(soj{i}>10^(-12))):end);
%     waitold{i}=deconv(stemp,prob);
% end    

% compute overall waiting time
waitT=zeros(size(wait{1}));
for i=1:K
   waitT(1,1:size(wait{i},2))=waitT(1,1:size(wait{i},2))+lambdas(i)*wait{i}; 
end
waitT=waitT/sum(lambdas);

% compute queue length per type
piTA=zeros(m*K,size(pi,2));
for k=1:K
    for j=1:m
        piTA(m*(k-1)+j,:)=sum(pi(m*smk(k)+j:m:m*smk(k+1),:),1);
        % elements piTA(m*(k-1)+1:m*k,i) give type k in service
    end    
end
piA=zeros(m,size(pi,2));
for j=1:m
    piA(j,:)=sum(piTA(j:m:end,:),1); 
end 
Dtemps=zeros(m,size(pi,2));
for k=1:K
    piAk=piA-piTA(m*(k-1)+1:m*k,:); % no type k in service
    ql{k}=zeros(1,size(pi,2)+1);
    ql{k}(1)=sum(pi0);
    Dnok=Dsum-D{k};
    Dtemps(:,1)=ones(m,1);
    for i=1:size(pi,2)
        % no type k in service
        ql{k}(1:i)=ql{k}(1:i)+piAk(:,i)'*Dtemps(:,1:i);
        % type k in service
        ql{k}(2:i+1)=ql{k}(2:i+1)+piTA(m*(k-1)+1:m*k,i)'*Dtemps(:,1:i);
        Dtemps(:,1:(i+1))=[Dnok*Dtemps(:,1:i) zeros(m,1)]+...
            [zeros(m,1) D{k}*Dtemps(:,1:i)];
    end
    ql{k}=ql{k}(1:max(find(ql{k}>10^(-10))));
end

% compute overall queue length
Dtemps=zeros(m,size(pi,2));
Dtemps(:,1)=ones(m,1);
qlT=zeros(size(pi,2)+1);
qlT(1)=sum(pi0);
D1=Dsum-D0;
for i=1:size(pi,2)
    % customer in service
    qlT(2:i+1)=qlT(2:i+1)+piA(:,i)'*Dtemps(:,1:i);
    Dtemps(:,1:(i+1))=[D0*Dtemps(:,1:i) zeros(m,1)]+...
        [zeros(m,1) D1*Dtemps(:,1:i)];
end
qlT=qlT(1:max(find(qlT>10^(-10))));