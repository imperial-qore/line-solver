function [ql,soj_alpha,wait_alpha,Smat]=Q_CT_MAP_MAP_1(C0,C1,D0,D1,varargin)
%[ql,soj_alpha,wait_alpha,S]=Q_CT_MAP_MAP_1(C0,C1,D0,D1) 
%   computes the  Sojourn time and Waiting time distribution of 
%   a continuous time MAP/MAP/1/FCFS queue 
%   
%   INPUT PARAMETERS:
%   * MAP arrival process (with m_a states)
%     the m_axm_a matrices C0 and C1 characterize the MAP arrival process
%
%   * MAP service process (with m_s states)
%     the m_sxm_s matrices D0 and D1 characterize the MAP service process
%
%   RETURN VALUES: 
%   * Queue length distribution   
%     ql(i) = Prob[(i-1) customers in the queue]
%   * Sojourn time distribution, 
%     is PH characterized by (soj_alpha,Smat)
%   * Waiting time distribution, 
%     is PH characterized by (wait_alpha,Smat)
%
%   USES: QBD Solver and QBD_pi from SMCSolver and Q_Sylvest
%
%   OPTIONAL PARAMETERS:
% 
%       Mode: Specifies the underlying function to compute the R matrix of 
%           the underlying QBD can be selected using the following 
%           parameter values (default: 'CR')
%               'CR' : Cyclic Reduction [Bini, Meini]
%               'FI' : Functional Iterations [Neuts]
%               'IS' : Invariant Subspace [Akar, Sohraby]
%               'LR' : Logaritmic Reduction [Latouche, Ramaswami]
%               'NI' : Newton Iteration
%       
%       MaxNumComp: Maximum number of components for the vectors containig
%       the performance measure.
%       
%       Verbose: When set to 1, the computation progress is printed
%       (default:0).
%       
%       Optfname: Optional parameters for the underlying function fname.
%       These parameters are included in a cell with one entry holding
%       the name of the parameter and the next entry the parameter
%       value. In this function, fname can be equal to:
%           'QBD_CR' : Options for Cyclic Reduction [Bini, Meini]
%           'QBD_FI' : Options for Functional Iterations [Neuts]
%           'QBD_IS' : Options for Invariant Subspace [Akar, Sohraby]
%           'QBD_LR' : Options for Logaritmic Reduction [Latouche, Ramaswami]
%           'QBD_NI' : Options for Newton Iteration


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


OptionValues{1}=['CR      ';
                 'FI      ';
                 'IS      ';
                 'LR      ';
                 'NI      '];
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
Q_CT_MAP_ParsePara(C0,'C0',C1,'C1')
Q_CT_MAP_ParsePara(D0,'D0',D1,'D1')

% Parse Optional Parameters
options=Q_ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);
% Check for unused parameter
Q_CheckUnusedParaQBD(options);


% Determine constants
ma=size(C0,1);
ms=size(D0,1);
mtot=ma*ms;

% Test the load of the queue
invC0=(-C0)^(-1);
pi_a=stat(C1*invC0);
lambda=sum(pi_a*C1);
invD0=(-D0)^(-1);
pi_s=stat(D1*invD0);
mu=sum(pi_s*D1);

load=lambda/mu;
if load >= 1
    error('MATLAB:Q_CT_MAP_MAP_1:LoadExceedsOne',...
                        'The load %d of the system exceeds one',load);
end  

%  Compute Queue length
% Compute classic QBD blocks A0, A1 and A2
Am1 = kron(eye(ma),D1);
A0 = kron(eye(ma),D0)+kron(C0,eye(ms));
A1 = kron(C1,eye(ms));
B0 = kron(C0,eye(ms));

if  (strfind(options.Mode,'FI')>0)
    [G,R,U]=QBD_FI(Am1,A0,A1,options.OptQBD_FI{:});
elseif (strfind(options.Mode,'LR')>0)
    [G,R,U]=QBD_LR(Am1,A0,A1,options.OptQBD_LR{:});
elseif (strfind(options.Mode,'IS')>0)
    [G,R,U]=QBD_IS(Am1,A0,A1,options.OptQBD_IS{:});
elseif (strfind(options.Mode,'NI')>0)
    [G,R,U]=QBD_NI(Am1,A0,A1,options.OptQBD_NI{:});
else
    [G,R,U]=QBD_CR(Am1,A0,A1,options.OptQBD_CR{:});
end

stv = QBD_pi(Am1,B0,R,'MaxNumComp',options.MaxNumComp,'Verbose',options.Verbose);

ql = zeros(1,size(stv,2)/(mtot));
for i=1:size(ql,2)
    ql(i) = sum(stv((i-1)*mtot+1:i*mtot));
end


% Compute Sojourn and Waiting time PH representation
if(nargout > 1)
    % Compute T directly
    Tnew=kron(eye(ma),D0)+(-U)^(-1)*kron(C1,D1);
    
    % Compute Smat
    theta_tot=kron(pi_a*C1,pi_s/mu)/load;
    nonz=find(theta_tot>0);
    theta_tot_red=theta_tot(nonz);
    Tnew=Tnew(:,nonz);
    Tnew=Tnew(nonz,:);
    Smat=diag(theta_tot_red)^(-1)*Tnew'*diag(theta_tot_red);

    % alpha vector of PH representation of Sojourn time
    soj_alpha=load*(diag(theta_tot)*kron(ones(ma,1),sum(D1,2)))'/lambda;
    soj_alpha=soj_alpha(nonz);

    % alpha vector of PH representation of Waiting time
    wait_alpha=load*(diag(theta_tot)*(-U)^(-1)*kron(sum(C1,2),sum(D1,2)))'/lambda;
    wait_alpha=wait_alpha(nonz);
end