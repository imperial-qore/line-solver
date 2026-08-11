function [wait,waitT,soj,sojT]=Q_DT_MMAPK_SMK_1(D0,D,C,varargin)
%   [wait,waitT,soj,sojT]=Q_DT_MMAPK_SMK_1(D0,D,C) computes the
%   Waiting Time and Sojourn time distribution of a 
%   Discrete-Time MMAP[K]/SM[K]/1/FCFS queue (per type and overall)
%   
%   INPUT PARAMETERS:
%   * MMAP[K] arrival process (with m_a states)
%     D{i} holds the m_axm_a matrix D_i, together with the mxm matrix D0 
%     these characterize the MMAP[K] arrival process
%
%   * SM[K] service time distributions (i=1:K) with m_s states
%     C{i} holds the matrix [C_i(1) C_i(2)...C_i(tmax(i))], where the 
%     m_sxm_s matrix C_i(j) holds the transition probabilities that the 
%     service time of a type-i customer takes j time units (j=1:tmax(i))
%
%   RETURN VALUES:
%   * Per type queue length distribution, 
%     ql{k}(i) = Prob[(i-1) type k customers in the queue]
%   * Overall queue length distribution,
%     qlT(i) = Prob[(i-1) customers in the queue (any type)]
%   * Per type waiting time distribution, 
%     wait{k}(i) = Prob[a type k customer has waiting time = (i-1)]
%   * Overall waiting time distribution,
%     waitT(i) = Prob[a customer (of any type) has waiting time = (i-1)]
%   * Per type sojourn time distribution, 
%     soj{k}(i) = Prob[a type k customer has sojourn time = (i-1)]
%   * Overall sojourn time distribution,
%     sojT(i) = Prob[a customer (of any type) has sojourn time = (i-1)]
%
%   OPTIONAL PARAMETERS:
%       Mode: The underlying function to compute the G matrix of the 
%             underlying MG1-type Markov chain can be selected using 
%             the following parameter values (default: 'CR')
%               
%               'CR' : Cyclic Reduction [Bini, Meini]
%               'FI' : Functional Iterations [Neuts]
%               'IS' : Invariant Subspace [Akar, Sohraby]
%               'RR' : Ramaswami Reduction [Bini,Meini,Ramaswami]
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
%               'MG1_CR' : Options for Cyclic Reduction [Bini, Meini]
%               'MG1_FI' : Options for Functional Iterations [Neuts]
%               'MG1_IS' : Options for Invariant Subspace [Akar, Sohraby]
%               'MG1_RR' : Options for ramaswami Reduction [Bini,Meini,Ramaswami]
%
%   USES: MG1 Solver and MG1_pi of the SMCSolver tool


OptionNames=[
             'Mode              ';   
             'MaxNumComp        ';
             'Verbose           ';
             'OptMG1_CR         '; 
             'OptMG1_FI         '; 
             'OptMG1_IS         '; 
             'OptMG1_RR         '];
         
OptionTypes=[
             'char   ';
             'numeric';
             'numeric';
             'cell   '; 
             'cell   '; 
             'cell   '; 
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
options.OptMG1_CR=cell(0);
options.OptMG1_FI=cell(0);
options.OptMG1_IS=cell(0);
options.OptMG1_RR=cell(0);

% Parse Parameters
K=size(D,2);
Q_DT_MMAPK_ParsePara(D0,'D0',D,'D');
Q_DT_SMK_ParsePara(C,'C',0);

% Parse Optional Parameters
options=Q_ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

% Parse Optional Parameters
Q_CheckUnusedParaMG1(options);

% Determine constants
ma=size(D0,1);
ms=size(C{1},1);
for k = 1:K
    tmax(k) = size(C{k},2)/ms;
end
mtot=ma*ms;

% Test the load of the queue
Dplus = D{1};
for i=2:K
    Dplus=Dplus+D{i};
end
Dsum=D0+Dplus;
thetaA=stat(Dsum);
Csum = reshape(sum(reshape(C{1}, ms*ms, tmax(1)), 2), ms, ms);
thetaS=stat(Csum);

lambdas=zeros(1,K);
means=zeros(1,K);
for k=1:K
    lambdas(k)=thetaA*D{k}*ones(ma,1);
    Cws = reshape(reshape(C{k},ms*ms, tmax(k))*[1:tmax(k)]', ms, ms);
    means(k)=thetaS*Cws*ones(ms,1);
end   
load=lambdas*means';
if load >= 1
    error('MATLAB:Q_DT_MMAPK_SMK_1:LoadExceedsOne',...
                        'The load %d of the system exceeds one',load);
end    

% Construct MG1 blocks
A0 = kron(D0,eye(ms));
A = [A0 zeros(mtot,max(tmax)*mtot)];
for k = 1:K
    A(:,mtot+1:2*mtot) = A(:,mtot+1:2*mtot) + kron( D{k},C{k}(:,1:ms) );
    for j = 2:tmax(k)
        A(:,j*mtot+1:(j+1)*mtot) = A(:,j*mtot+1:(j+1)*mtot) + kron( D{k},C{k}(:,(j-1)*ms+1:j*ms) );
    end
end

%SMC Solver
switch options.Mode
    case 'CR'
        G = MG1_CR(A,options.OptMG1_CR{:});
    case 'FI'
        G = MG1_FI(A,options.OptMG1_FI{:});
    case 'RR'
        G = MG1_RR(A,options.OptMG1_RR{:});
    case 'IS'
        G = MG1_IS(A,options.OptMG1_IS{:});
end
pi=MG1_pi([],A,G,'MaxNumComp',options.MaxNumComp,'Verbose', options.Verbose);


% compute overall waiting time distribution
lambdaA = sum(lambdas);
Dplus_vec = kron(sum(Dplus,2),ones(ms,1));
D_vec = cell(1,K);
for i=1:K
    D_vec{i} = kron(sum(D{i},2),ones(ms,1));
end
sizeW = size(pi,2)/mtot-1;
waitT=zeros(1,sizeW);
waitT(1) = (pi(1:mtot)+pi(mtot+1:2*mtot))*Dplus_vec/lambdaA;
for i = 2:sizeW
    waitT(i) = pi(i*mtot+1:(i+1)*mtot)*Dplus_vec/lambdaA;
end
waitT = waitT(1:max(find(waitT>10^(-10))));

% compute per-type waiting time distribution
wait = cell(1,K);
for k = 1:K
    waitTemp=zeros(1,sizeW);
    waitTemp(1) = (pi(1:mtot)+pi(mtot+1:2*mtot))*D_vec{k}/lambdas(k);
    for i = 2:sizeW
        waitTemp(i) = pi(i*mtot+1:(i+1)*mtot)*D_vec{k}/lambdas(k);
    end
    wait{k} = waitTemp;
    wait{k} = wait{k}(1:max(find(wait{k}>10^(-10))));
end    

% compute overall sojourn time distribution
sizeS = size(pi,2)/mtot + max(tmax);
sojT = zeros(1,sizeS);
Ctilde = zeros(mtot,max(tmax)); %Ctilde(n) = (Dplus kron C(n)) ones(mtot,1)
Ctilde_k=cell(1,K);
for k = 1:K
    Ctemp = reshape(sum(reshape(C{k}', ms,ms*tmax(k)),1), tmax(k), ms)';
    Ctemp = kron( sum(D{k},2), Ctemp );
    Ctilde_k{k} = Ctemp;
    Ctilde(:,1:tmax(k)) = Ctilde(:,1:tmax(k)) + Ctilde_k{k};
end




for n = 1:sizeS
    sojTemp = 0;
    if n <= max(tmax)
        sojTemp = pi(1:mtot)*Ctilde(:,n);
    end
    sojTemp = sojTemp + pi(max(1,n+1-max(tmax))*mtot+1:(min(n,size(pi,2)/mtot-1)+1)*mtot)...
        *reshape(Ctilde(:,n+1-[max(1,n+1-max(tmax)):min(n,size(pi,2)/mtot-1)]),[],1);
    sojT(n) = sojTemp/lambdaA;
end
sojT = [0 sojT];
sojT = sojT(1:max(find(sojT>10^(-10))));

% compute per-type sojourn time distribution
soj = cell(1,K);
for k = 1:K
    sizeS_k = size(pi,2)/mtot + tmax(k);
    soj{k} = zeros(1,sizeS_k);
    for n = 1:sizeS
        sojTemp = 0;
        if n <= tmax(k)
            sojTemp = pi(1:mtot)*Ctilde_k{k}(:,n);
        end
        sojTemp = sojTemp + pi(max(1,n+1-tmax(k))*mtot+1:(min(n,size(pi,2)/mtot-1)+1)*mtot)...
            *reshape(Ctilde_k{k}(:,n+1-[max(1,n+1-tmax(k)):min(n,size(pi,2)/mtot-1)]),[],1);
        soj{k}(n) = sojTemp/lambdas(k);
    end
    soj{k} = [0 soj{k}];
    soj{k}=soj{k}(1:max(find(soj{k}>10^(-10))));
end