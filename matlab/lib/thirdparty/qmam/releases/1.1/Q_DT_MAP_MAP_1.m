function [ql,wait,soj]=Q_DT_MAP_MAP_1(C0,C1,D0,D1,varargin)
%   [ql,wait,soj] = Q_DT_MAP_MAP_1(C0,C1,D0,D1) computes the Queue length,
%   Sojourn time and Waiting time distribution of a discrete-time 
%   MAP/MAP/1/FCFS queue
%   
%   INPUT PARAMETERS:
%   * MAP arrival process (with m_a states)
%     the m_axm_a matrices C0 and C1 characterize MAP arrival process
%
%   * MAP service process (with m_s states)
%     the m_sxm_s matrices D0 and D1 characterize MAP service process
%
%   RETURN VALUES:
%   * Queue length distribution, 
%     ql(i) = Prob[(i-1) customers in the queue]
%   * Waiting time distribution,
%     wait(i) = Prob[a customer has waiting time = (i-1)]
%   * Sojourn time distribution,
%     soj(i) = Prob[a customer has sojourn time = (i-1)]
%
%   OPTIONAL PARAMETERS:
%       Mode: The underlying function to compute the R matrix of the 
%           underlying QBD can be selected using the following 
%           parameter values (default: 'CR')
%               'CR' : Cyclic Reduction [Bini, Meini]
%               'FI' : Functional Iterations [Neuts]
%               'IS' : Invariant Subspace [Akar, Sohraby]
%               'LR' : Logaritmic Reduction [Latouche, Ramaswami]
%               'NI' : Newton Iteration
%
%       MaxNumComp: Maximum number of components for the vectors containig
%           the performance measure.
%       
%       Verbose: When set to 1, the computation progress is printed
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
Q_DT_MAP_ParsePara(C0,'C0',C1,'C1');
Q_DT_MAP_ParsePara(D0,'D0',D1,'D1');

% Parse Optional Parameters
options=Q_ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

% Parse Optional Parameters
Q_CheckUnusedParaQBD(options);


% arrivals
ma = size(C0,1);
C = C0+C1;
avga = stat(C)*C1*ones(ma,1);

% service
ms = size(D0,1);
D = D0+D1;
piD = stat(D);
avgs = piD*D1*ones(ms,1);

mtot = ms*ma;
rho = avga/avgs;
if rho >= 1
    error('MATLAB:Q_DT_MAP_MAP_1:LoadExceedsOne',...
                        'The load %d of the system exceeds one',rho);
end    

% Compute classic QBD blocks A0, A1 and A2
Am1 = kron(C0,D1);
A0 = kron(C0,D0)+kron(C1,D1);
A1 = kron(C1,D0);
Bm1 = kron(C0,D1);
B0 = kron(C0,eye(ms));
B1 = kron(C1,eye(ms));

if  (strfind(options.Mode,'FI')>0)
    [G,R]=QBD_FI(Am1,A0,A1,options.OptQBD_FI{:});
elseif (strfind(options.Mode,'LR')>0)
    [G,R]=QBD_LR(Am1,A0,A1,options.OptQBD_LR{:});
elseif (strfind(options.Mode,'IS')>0)
    [G,R]=QBD_IS(Am1,A0,A1,options.OptQBD_IS{:});
elseif (strfind(options.Mode,'NI')>0)
    [G,R]=QBD_NI(Am1,A0,A1,options.OptQBD_NI{:});
else
    [G,R]=QBD_CR(Am1,A0,A1,options.OptQBD_CR{:});
end

stv = QBD_pi(Bm1,B0,R,'Boundary',[B1; A0+R*Am1],'MaxNumComp',options.MaxNumComp,'Verbose',options.Verbose);

% Compute queue length distribution
ql = zeros(1,size(stv,2)/mtot);
for i=1:size(ql,2)
    ql(i) = sum(stv((i-1)*mtot+1:i*mtot));
end

% Compute Waiting & Sojourn time distribution
if(nargout > 1)
    % 1. (number in queue, service phase) after arrival
    stva=reshape(stv,mtot,size(stv,2)/mtot)';
    stva=stva(2:end,:)*kron(sum(C1,2),D1)+...
        [stva(1,:)*kron(sum(C1,2),eye(ms)); ...
         stva(2:end-1,:)*kron(sum(C1,2),D0)];
    stva=stva/sum(sum(stva));
    len=size(stva,1);
    stva=reshape(stva',1,len*ms);
    % 2. some memory preallocation
    wait=zeros(1,len*2*ceil(1/avgs));
    soj=zeros(1,len*2*ceil(1/avgs));

    temp=zeros(ms*len,1);
    temp(1:ms)=sum(D1,2);
    wait(1)=sum(stva(1:ms),2);  % waiting time = 0
    soj(1)=0;  % sojourn time = 0
    i=2;j=2;
    while (min(sum(wait),sum(soj))<1-10^(-10))
        wait(j)=stva(ms+1:i*ms)*temp(1:(i-1)*ms);
        soj(j)=stva(1:(i-1)*ms)*temp(1:(i-1)*ms);
        temp=reshape(temp,ms,len);
        temp(:,1:i)=D0*[temp(:,1:i-1) zeros(ms,1)]+D1*[zeros(ms,1) temp(:,1:i-1)];
        temp=reshape(temp,ms*len,1);
        i=min(i+1,len);
        j=j+1;
    end    
    wait=wait(1:max(find(wait>10^(-10))));
    soj=soj(1:max(find(soj>10^(-10))));
end