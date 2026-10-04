function [ql,wait,soj]=Q_DT_PH_PH_1(alpha,T,beta,S,varargin)
%[ql,wt]=Q_DT_PH_PH_1(alpha,T,beta,S) computes the Queuelength, 
%   Waiting time and Sojourn time distribution of a 
%   Discrete-Time PH/PH/1/FCFS queue
%   
%   INPUT PARAMETERS:
%   * PH inter-arrival time distribution
%     alpha, the 1xma vector of the PH arrival process
%     T, the maxma matrix of the PH arrival process
%
%   * PH service time distributions
%     beta, the 1xms vector of the PH service time 
%     S, the msxms matrix of the PH service time
%
%   RETURN VALUES:
%   * Queue length distribution, 
%     ql(i) = Prob[(i-1) customers in the queue]
%   * Waiting time distribution, 
%     wt(i) = Prob[a customer has waiting time = (i-1)]
%   * Sojourn time distribution,
%     soj(i) = Prob[a customer has sojourn time = (i-1)]
%
%   OPTIONAL PARAMETERS:
%       MaxNumComp: Maximum number of components for the vectors containig
%           the performance measure.
%       
%       Verbose: When set to 1, the computation progress is printed
%           (default:0).


OptionNames=[
             'MaxNumComp        ';
             'Verbose           '];
OptionTypes=[
             'numeric';
             'numeric'];

OptionValues=cell(0);
 
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.MaxNumComp = 1000;
options.Verbose = 0;

% Parse Parameters
Q_DT_PH_ParsePara(alpha,'alpha',T,'T');
Q_DT_PH_ParsePara(beta,'beta',S,'S');

% Parse Optional Parameters
options=Q_ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

% arrivals
ma = size(alpha,2);
t = ones(ma,1)-T*ones(ma,1);
avga = (alpha*inv(eye(ma)-T)*ones(ma,1));

% service
ms = size(beta,2);
s = ones(ms,1)-S*ones(ms,1);
avgs = beta*inv(eye(ms)-S)*ones(ms,1);

mtot = ma*ms;
rho = avgs/avga;

if rho >= 1
    error('MATLAB:Q_DT_PH_PH_1:LoadExceedsOne',...
                        'The load %d of the system exceeds one',rho);
end    

% Compute classic QBD blocks B0, B1, A0, A1 and A2
A0 = kron(t*alpha,S);
A1 = kron(T,S)+kron(t*alpha,s*beta);
% A2 = kron(T,s*beta);
B0 = kron(t*alpha,eye(ms));
% B1 = kron(T,eye(ms));

% Compute QBD blocks in approach Latouche & Ramaswami
A0pp = kron(alpha,eye(ms))*inv(eye(mtot)-kron(T,S))*kron(t,S);
A0mp = kron(eye(ma),beta)*inv(eye(mtot)-kron(T,S))*kron(t,S);
A2pm = kron(alpha,eye(ms))*inv(eye(mtot)-kron(T,S))*kron(T,s);
A2mm = kron(eye(ma),beta)*inv(eye(mtot)-kron(T,S))*kron(T,s);
A1pp = kron(alpha,eye(ms))*inv(eye(mtot)-kron(T,S))*kron(t,s*beta);
A1mp = kron(eye(ma),beta)*inv(eye(mtot)-kron(T,S))*kron(t,s*beta);

A0n = [A0pp         zeros(ms,ma);   A0mp            zeros(ma)];
A2n = [zeros(ms)    A2pm;           zeros(ma,ms)    A2mm];
A1n = [A1pp         zeros(ms,ma);   A1mp            zeros(ma)];


% Compute Gamma
itB0 = inv(eye(ma+ms)-A1n)*A0n;
itB2 = inv(eye(ma+ms)-A1n)*A2n;
Gamma = itB2(1:ms,ms+1:end);
itT = itB0;
while norm(ones(ms,1)-Gamma*ones(ma,1)) > 1e-10
        itA1 = itB0*itB2 + itB2*itB0;
        itB0 = inv(eye(ma+ms)-itA1)*itB0^2;
        itB2 = inv(eye(ma+ms)-itA1)*itB2^2;
        tmp = itT*itB2;
        Gamma = Gamma + tmp(1:ms,ms+1:end);
        itT = itT*itB0;
end
Gm = inv(eye(ma)-A0mp*Gamma)*A2mm + inv(eye(ma)-A0mp*Gamma)*A1mp*Gamma;
 
% Compute queue length distribution
Gstar = inv(eye(mtot)-A1)*(kron(T,s*beta) + kron(t,S)*Gamma*Gm*kron(eye(ma),beta));
Rstar = A0*inv(eye(mtot)-A1-A0*Gstar);
Rstar0 = B0*inv(eye(mtot)-A1-A0*Gstar);

% Compute pi_0
stv = kron((1-rho)*inv(beta*Gamma*inv(eye(ma)-T)*ones(ma,1))*beta*Gamma*inv(eye(ma)-T),beta);
% Compute pi_1
stv(2,1:mtot) = stv*Rstar0;
% Compute pi_2,...
sumpi=sum(sum(stv));
numit=2;
while (sumpi < 1-10^(-10) && numit < 1+options.MaxNumComp)
    stv(numit+1,1:mtot)=stv(numit,:)*Rstar; %compute pi_(numit+1)
    numit=numit+1;
    sumpi=sumpi+sum(stv(numit,:));
    if (~mod(numit,options.Verbose))
        fprintf('Accumulated mass after %d iterations: %d\n',numit,sumpi);
        drawnow;
    end
end  
ql = sum(stv,2)';
stv=reshape(stv',1,[]); 
if (numit == 1+options.MaxNumComp)
    warning('Maximum Number of Components %d reached',numit-1);
end


% Compute Waiting & Sojourn time distribution
if(nargout > 1)
    % 1. (number in queue, service phase) after arrival
    stva=reshape(stv,mtot,size(stv,2)/mtot)';
    stva=stva(2:end,:)*kron(t,s*beta)+...
        [stva(1,:)*kron(t,eye(ms)); ...
         stva(2:end-1,:)*kron(t,S)];
    stva=stva/sum(sum(stva));
    len=size(stva,1);
    stva=reshape(stva',1,len*ms);
    % 2. some memory preallocation
    wait=zeros(1,len*2*ceil(1/avgs));
    soj=zeros(1,len*2*ceil(1/avgs));

    temp=zeros(ms*len,1);
    temp(1:ms)=s;
    wait(1)=sum(stva(1:ms),2);  % waiting time = 0
    soj(1)=0;  % sojourn time = 0
    i=2;j=2;
    while (min(sum(wait),sum(soj))<1-10^(-10))
        wait(j)=stva(ms+1:i*ms)*temp(1:(i-1)*ms);
        soj(j)=stva(1:(i-1)*ms)*temp(1:(i-1)*ms);
        temp=reshape(temp,ms,len);
        temp(:,1:i)=S*[temp(:,1:i-1) zeros(ms,1)]+s*beta*[zeros(ms,1) temp(:,1:i-1)];
        temp=reshape(temp,ms*len,1);
        i=min(i+1,len);
        j=j+1;
    end    
    wait=wait(1:max(find(wait>10^(-10))));
    soj=soj(1:max(find(soj>10^(-10))));
end