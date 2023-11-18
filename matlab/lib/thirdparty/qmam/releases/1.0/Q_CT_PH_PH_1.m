function [ql,wait_alpha,wait_T]=Q_CT_PH_PH_1(alpha,T,beta,S,varargin)
%[ql,wtalpha,wtT]=Q_CT_PH_PH_1(alpha,T,beta,S) computes the Queuelength 
%   and Waiting time distribution of a Continuous-Time PH/PH/1/FCFS queue
%   
%   INPUT PARAMETERS:
%   * PH inter-arrival time distribution
%     alpha: the 1xma vector of the PH arrival process
%     T: the maxma matrix of the PH arrival process
%
%   * PH service time distributions
%     beta: the 1xms vector of the PH service time 
%     S: the msxms matrix of the PH service time
%
%   RETURN VALUES:
%   * Queue length distribution, 
%     ql(i) = Prob[(i-1) customers in the queue]
%   * Waiting time distribution, 
%     is PH characterized by (wait_alpha, wait_T)
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
Q_CT_PH_ParsePara(alpha,'alpha',T,'T');
Q_CT_PH_ParsePara(beta,'beta',S,'S');

% Parse Optional Parameters
options=Q_ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);


% arrivals
ma = size(alpha,2);
t = -T*ones(ma,1);
avgt = (alpha*inv(-T)*ones(ma,1));

% service
ms = size(beta,2);
s = -S*ones(ms,1);
avgs = beta*inv(-S)*ones(ms,1);

mtot = ms*ma;
rho = avgs/avgt;
if rho >= 1
    error('MATLAB:Q_CT_PH_PH_1:LoadExceedsOne',...
                        'The load %d of the system exceeds one',rho);
end    

% Compute classic QBD blocks A0, A1 and A2
A0 = kron(t*alpha,eye(ms));
A1 = kron(T,eye(ms))+kron(eye(ma),S);
% A2 = kron(eye(ma),s*beta);
% B = kron(T,eye(ms));

% Compute QBD blocks in approach Latouche & Ramaswami
invmA1 = inv(-A1);
A0pp = kron(alpha,eye(ms))*invmA1*kron(t,eye(ms));
A0mp = kron(eye(ma),beta)*invmA1*kron(t,eye(ms));
A2pm = kron(alpha,eye(ms))*invmA1*kron(eye(ma),s);
A2mm = kron(eye(ma),beta)*invmA1*kron(eye(ma),s);

A0n = [A0pp zeros(ms,ma); A0mp zeros(ma)];
A2n = [zeros(ms) A2pm;zeros(ma,ms) A2mm];

% Compute matrix Gamma: NE corner of matrix G
itB0 = A0n;
itB2 = A2n;
Gamma = itB2(1:ms,ms+1:end);
itT = itB0;
check=1;
numit = 1;
while check > 10e-14
    itA1 = itB0*itB2 + itB2*itB0;
    itB0 = inv(eye(ma+ms)-itA1)*itB0^2;
    itB2 = inv(eye(ma+ms)-itA1)*itB2^2;
    tmp = itT*itB2;
    Gamma = Gamma + tmp(1:ms,ms+1:end);
    itT = itT*itB0;
    check = norm(ones(ms,1)-Gamma*ones(ma,1));
    numit=numit+1;
    if (options.Verbose==1)
        fprintf('Check after %d iterations: %d\n',numit,check);
        drawnow;
    end
end

Gm = inv(eye(ma)-A0mp*Gamma)*A2mm;
R_Gam = A0pp*inv(eye(ms)-Gamma*A0mp);


% Compute queue length distribution
Gstar = invmA1*kron(eye(ma),s*beta) + invmA1*kron(t,eye(ms))*Gamma*Gm*kron(eye(ma),beta);
Rstar = A0*inv(-A1-A0*Gstar);

% Compute pi_0
pi = kron((1-rho)*inv(beta*Gamma*inv(-T)*ones(ma,1))*beta*Gamma*inv(-T),beta);
% Compute pi_1,...
sumpi=sum(pi);
numit=1;
while (sumpi < 1-10^(-10) && numit < 1+options.MaxNumComp)
    pi(numit+1,1:mtot)=pi(numit,:)*Rstar; %compute pi_(numit+1)
    numit=numit+1;
    sumpi=sumpi+sum(pi(numit,:));
    if (~mod(numit,options.Verbose))
        fprintf('Accumulated mass after %d iterations: %d\n',numit,sumpi);
        drawnow;
    end
end
ql = sum(pi,2)';
pi=reshape(pi',1,[]);
if (numit == 1+options.MaxNumComp)
    warning('Maximum Number of Components %d reached',numit-1);
end


% Compute waiting time PH representation
if nargout > 1
    sigtilde = inv(beta*inv(S)*ones(ms,1))*beta*inv(S);
    Delta = diag(sigtilde);
    wait_T = inv(Delta)*(S+R_Gam*s*beta)'*Delta;
    sigrho = rho*sigtilde;
    theta = (-beta*inv(S)*ones(ms,1))*s'*Delta;
    D = inv(Delta)*R_Gam'*Delta;
    wait_alpha = (1-inv(beta*inv(eye(ms)-R_Gam)*ones(ms,1))*beta*ones(ms,1))*inv(theta*D*ones(ms,1))*theta*D;
end