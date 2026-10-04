function [ql,wait_alpha,Smat]=Q_CT_MAP_M_C(D0,D1,mu,c,varargin)
%   [ql,wait_alpha,Smat]=Q_CT_MAP_M_C(D0,D1,mu,c) computes the Queuelength 
%   and waiting time distribution of a MAP/M/c/FCFS queue 
%   
%   INPUT PARAMETERS:
%   * MAP(D0,D1) arrival process (with m states)
%
%   * mu exponential service rate
%     
%   * c number of servers
%
%   RETURN VALUES:
%   * Queue length distribution, 
%     ql(i) = Prob[(i-1) customers in the queue]
%
%   * Waiting time distribution, 
%     is PH characterized by (wait_alpha,Smat)
%
%   USES: USES: QBD Solver from SMCSolver and Q_Sylvest
%
%   OPTIONAL PARAMETERS:
% 
%       Mode: The underlying methods used to compute the performance
%           measures. For this function two different modes can be 
%           specified. The parameter value may include a specific option 
%           for one of them, or a combination of one name for each option, 
%           i.e., both 'Sylves', 'FI' and 'SylvesFI' are possible values 
%           for this parameter.
%
%           1. The method to solve the linear system at each iteration of
%           the algorithm to compute matrix T (default: 'Sylves')
%               'Sylves' : solves a Sylvester matrix equation at each step 
%                          using a Hessenberg algorithm
%               'Direct' : solves the Sylvester matrix equation at each
%                          step by rewriting it as a (large) system of 
%                          linear equations.
%
%           2. The underlying function to compute the R matrix of the 
%           underlying QBD can be selected using the following 
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
%       Verbose: When set to 1, the progress of the computation is printed
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


OptionValues{1}=['Direct  '; 
                 'Sylves  ';
                 'CR      ';
                 'FI      ';
                 'IS      ';
                 'LR      ';
                 'NI      ';
                 'DirectCR';
                 'DirectFI'; 
                 'DirectIS'; 
                 'DirectLR'; 
                 'DirectNI';
                 'SylvesCR';
                 'SylvesFI';
                 'SylvesIS';
                 'SylvesLR';
                 'SylvesNI'];
  
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.Mode='SylvesCR';
options.MaxNumComp = 1000;
options.Verbose = 0;
options.OptQBD_CR=cell(0);
options.OptQBD_FI=cell(0);
options.OptQBD_IS=cell(0);
options.OptQBD_LR=cell(0);
options.OptQBD_NI=cell(0);

% Parse Parameters
if (~isnumeric(mu))
    error('MATLAB:Q_CT_MAP_M_C_ParsePara:InvalidInput',...
        'the service rate mu has to be numeric');
end    
if (~isnumeric(c))
    error('MATLAB:Q_CT_MAP_M_C_ParsePara:InvalidInput',...
        'the number of servers c has to be numeric');
end    
% check real
if (~isreal(mu))
    error('MATLAB:Q_CT_MAP_M_C_ParsePara:InvalidInput',...
        'the service rate mu has to be real');
end    
if (~isreal(c))
    error('MATLAB:Q_CT_MAP_M_C_ParsePara:InvalidInput',...
        'the number of servers c has to be real');
end 
if (ceil(c)-floor(c)>0)
    error('MATLAB:Q_CT_MAP_M_C_ParsePara:InvalidInput',...
        'the number of servers c has to be an integer');
end 
%check positivity
if (mu < 10E-14)
    error('MATLAB:Q_CT_MAP_M_C_ParsePara:InvalidInput',...
        'the service rate mu has to be strictly positive');
end 
if (c < 10E-14)
    error('MATLAB:Q_CT_MAP_M_C_ParsePara:InvalidInput',...
        'the number of servers c has to be strictly positive');
end 
%arrival process
Q_CT_MAP_ParsePara(D0,'D0',D1,'D1')


% Parse Optional Parameters
options=Q_ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);
% Check for unused parameter
Q_CheckUnusedParaQBD(options);


% Test the load of the queue
invD0=(-D0)^(-1);
pi_a=stat(D1*invD0);
lambda=sum(pi_a*D1);

load=lambda/(mu*c);
if load >= 1
    error('MATLAB:Q_CT_MAP_M_C:LoadExceedsOne',...
                        'The load %d of the system exceeds one',load);
end  

% Determine constants
m = size(D0,2);
A0 = c*mu*eye(m);
A1 = D0 - c*mu*eye(m);
A2 = D1;

%Compute Matrix R
if  (strfind(options.Mode,'FI')>0)
    [G,R]=QBD_FI(A0,A1,A2,options.OptQBD_FI{:});
elseif (strfind(options.Mode,'LR')>0)
    [G,R]=QBD_LR(A0,A1,A2,options.OptQBD_LR{:});
elseif (strfind(options.Mode,'IS')>0)
    [G,R]=QBD_IS(A0,A1,A2,options.OptQBD_IS{:});
elseif (strfind(options.Mode,'NI')>0)
    [G,R]=QBD_NI(A0,A1,A2,options.OptQBD_NI{:});
else
    [G,R]=QBD_CR(A0,A1,A2,options.OptQBD_CR{:});
end

%Gaver-Jacobs-Latouche LD-QBD
piGJL = zeros(1,c*m);
if c > 1
    invC = cell(1,c-1);
    invC{1} = inv(-D0);
    for i =2:c-1
        invC{i} = inv(-D0+(i-1)*mu*eye(m) - (i-1)*mu*invC{i-1}*D1);
    end

    piGJL((c-1)*m+1:c*m)  = stat(D0-(c-1)*mu*eye(m) + R*A0 + (c-1)*mu*invC{c-1}*D1 + eye(m));
    for i =c-1:-1:1
        piGJL((i-1)*m+1:i*m) = piGJL(i*m+1:(i+1)*m)*i*mu*invC{i};
    end
else
    piGJL(1:m)  = stat(D0 + R*A0 + eye(m));
end
K = sum(piGJL(1:(c-1)*m)) + piGJL((c-1)*m+1:c*m)*sum(inv(eye(m)-R),2);
piGJL = piGJL/K;

piC1 = piGJL((c-1)*m+1:c*m);
sumpi = sum(piGJL);


numit=1;
while (sumpi < 1-10^(-10) && numit < 1+options.MaxNumComp-c)
    piC1(numit+1,1:m)=piC1(numit,:)*R; 
    numit=numit+1;
    sumpi=sumpi+sum(piC1(numit,:));
    if (~mod(numit,options.Verbose))
        fprintf('Accumulated mass after %d iterations: %d\n',numit,sumpi);
        drawnow;
    end
end
qlC = sum(piC1,2);
piC1=reshape(piC1',1,[]);
piT = [piGJL(1:(c-1)*m) piC1];

% Queue length distribution
ql = zeros(1,c-1);
for i = 1:c-1
    ql(i) = sum(piGJL((i-1)*m+1:i*m));
end
ql = [ql qlC'];


% Waiting time distribution
if nargout > 1
    % a) Probability waiting = 0 and alpha vector
    prob_zero=sum(reshape(piGJL(1:c*m),m,c)'*sum(D1,2));
    prob_zero=prob_zero/sum(reshape(piT,m,size(piT,2)/m)'*sum(D1,2));
    temp=piGJL((c-1)*m+1:c*m)*(eye(m)-R)^(-1)*D1;
    alpha_vec=temp/sum(temp);

    % b) Ph representation
    % Compute Sojourn and Waiting time 
    Told=zeros(m,m);
    Tnew=-A0; % A0=eye(m)*mu*c
    if (strfind(options.Mode,'Direct')>0)
        while(max(max(abs(Told-Tnew)))>10^(-10))
            Told=Tnew;
            L=-reshape(eye(m),1,m^2)*(kron(Tnew',eye(m))+...
                kron(eye(m),D0))^(-1);
            Tnew=-A0+reshape(L,m,m)'*D1*mu*c;
        end
        L=reshape(L,m,m)';
    else
        [U,Tr]=schur(D0,'complex');
        %U'*Tr*U = A
        while(max(max(abs(Told-Tnew)))>10^(-10))
            Told=Tnew;
            L=Q_Sylvest(U,Tr,Tnew);
            Tnew=-A0+L*D1*mu*c;
        end
    end

    % Compute Smat
    theta_tot=alpha_vec;
    nonz=find(theta_tot>0);
    theta_tot_red=theta_tot(nonz);
    Tnew=Tnew(:,nonz);
    Tnew=Tnew(nonz,:);
    Smat=diag(theta_tot_red)^(-1)*Tnew'*diag(theta_tot_red);


    % alpha vector of PH representation of Waiting time
    rho_vec=sum(Tnew+A0,2);
    wait_alpha=alpha_vec.*rho_vec'/(alpha_vec*rho_vec);
    wait_alpha=(1-prob_zero)*wait_alpha(nonz);
end
