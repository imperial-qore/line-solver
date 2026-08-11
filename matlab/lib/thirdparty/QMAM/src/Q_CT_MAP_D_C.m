function [ql,w]=Q_CT_MAP_D_C(D0,D1,s,c,varargin)
%[ql]=Q_MAP_D_C(D0,D1,s,c) computes the Queuelength distribution 
%   of a MAP/D/c/FCFS queue with deterministic services of length s    
%   
%   INPUT PARAMETERS:
%   * MAP arrival process (with m states)
%     the mxm matrices D0 and D1 characterize the MAP arrival process
%     
%   * s length of deterministic service 
%
%   * c number of servers
%     
%   RETURN VALUES:
%   * Queue length distribution, 
%     ql(i) = Prob[(i-1) customers in the queue]
%   * Waiting time distribution, 
%     wt(i) = Prob[a customer has waiting time <= (i-1)*s/NumSteps]
%
%   USES: NSF_GHT and NSF_pi from the SMCSolver tool
%
%   OPTIONAL PARAMETERS:
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
%           'NSF_GHT': Options for Non-Skip-Free Markov Chains [Gail,
%                       Hantler, Taylor]
%
%       NumSteps: Number of points in the intervals {[(k-1)*s, k*s), k>0}
%       in which the waiting time distribution is evaluated.
%       (default:1).
    

OptionNames=[
             'MaxNumComp        ';
             'Verbose           ';
             'OptNSF_GHT        ';
             'NumSteps          '];
OptionTypes=[
             'numeric';
             'numeric';
             'cell   ';
             'numeric'];

OptionValues=cell(0);
 
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.MaxNumComp=1000;
options.Verbose=0;
options.OptNSF_GHT=cell(0);
options.NumSteps=1;

% Parse Parameters
if (~isnumeric(s))
    error('MATLAB:Q_CT_MAP_D_C_ParsePara:InvalidInput',...
        'the service time s has to be numeric');
end    
if (~isnumeric(c))
    error('MATLAB:Q_CT_MAP_D_C_ParsePara:InvalidInput',...
        'the number of servers c has to be numeric');
end    
% check real
if (~isreal(s))
    error('MATLAB:Q_CT_MAP_D_C_ParsePara:InvalidInput',...
        'the service time s has to be real');
end    
if (~isreal(c))
    error('MATLAB:Q_CT_MAP_D_C_ParsePara:InvalidInput',...
        'the number of servers c has to be real');
end 
if (ceil(c)-floor(c)>0)
    error('MATLAB:Q_CT_MAP_D_C_ParsePara:InvalidInput',...
        'the number of servers c has to be an integer');
end 
%check positivity
if (s < 10E-14)
    error('MATLAB:Q_CT_MAP_D_C_ParsePara:InvalidInput',...
        'the service time s has to be strictly positive');
end 
if (c < 10E-14)
    error('MATLAB:Q_CT_MAP_D_C_ParsePara:InvalidInput',...
        'the number of servers c has to be strictly positive');
end 
%arrival process
Q_CT_MAP_ParsePara(D0,'D0',D1,'D1')

% Parse Optional Parameters
options=Q_ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

% Determine constants
m = size(D0,2);
lambda = max(abs(diag(D0)));
P0 = D0/lambda + eye(m);
P1 = D1/lambda;

% Test the load of the queue
thetaA = stat(P0+P1);
lambdaA = thetaA*sum(D1,2);
epsilon = 10^-12;

load=lambdaA*s/c;
if load >= 1-epsilon
    error('MATLAB:Q_CT_MAP_D_C:LoadExceedsOne',...
                        'The load %d of the system exceeds one',load);
end    


% Compute NSF blocks: [P(0,s) P(1,s) ... P(max,s)]
%Zero arrivals
P0s = expm(D0*s);
Ptot = P0s;


%Determine number of Poisson terms (Pterms) to compute P(k,a)
j = 1;
h(j) = exp(-lambda*s)*lambda*s;
sumh = h(j) +  exp(-lambda*s);
while sumh < 1 - epsilon
    h = [h; h(j)*lambda*s/(j+1)];
    j = j + 1;
    sumh = sumh + h(j);
end
Pterms = j;

Kold = cell(1,Pterms);
Kold{1} = eye(m);
for j = 2:Pterms
    Kold{j} = Kold{j-1}*P0;
end

k = 1;
probCum = min(sum(Ptot,2));
while probCum < 1 - epsilon
    K{1} = P1*Kold{1};
    htemp = h(k);
    Ps{k} = htemp*K{1};
    for j = 2:Pterms
        K{j} = P0*K{j-1} + P1*Kold{j};
        htemp = htemp*lambda*s/(k-1+j);
        Ps{k} = Ps{k} + htemp*K{j};
    end
    Ptot = Ptot + Ps{k};
    
    for j = 1:Pterms
        Kold{j} = K{j};
    end
    probCum = min(sum(Ptot,2));
    k = k + 1;
end

%Compute matrix G
A = zeros(m,max(c+1,k)*m);
A(1:m,1:m) = P0s;
for i = 1:k-1
    A(1:m,i*m+1:(i+1)*m) = Ps{i};
end

G = NSF_GHT(A, c, options.OptNSF_GHT{:});
pi=NSF_pi([],A,G,'MaxNumComp', options.MaxNumComp, 'Verbose', options.Verbose);

% Compute queue length 
ql = zeros(1, size(pi,2)/m);
for i = 1:size(pi,2)/m
    ql(i)=sum(pi((i-1)*m+1:i*m));
end

% Waiting time distribution
if nargout > 1
    %Compute Waiting Time distribution W(t) for t = {0,s,2s,...}
    w(1) = sum(pi(1:m*c));
    WTacum = w(1);
    i = 2;
    while WTacum < 1 - 10^-10 && i*m*c < size(pi,2)
        w(i) = w(i-1) + sum( pi((i-1)*m*c+1:i*m*c) );
        WTacum = w(i);
        i = i+1;
    end
        
    
    %Compute Waiting Time distribution W(t) for t = {tau,s+tau,2s+tau,...}
    %with tau = i*s/NumSteps, for i ={1,2,...,NumSteps-1}
    options.NumSteps = floor(options.NumSteps)
    if options.NumSteps > 1
        for step = 1:options.NumSteps-1
            %General tau = i*s/NumSteps
            tau = step*s/options.NumSteps;
            
            %Compute P(k,tau), k\geq0
            %Determine number of Poisson terms (Pterms) to compute P(k,tau)
            clear h_tau K Kold Ptau;
            j = 1;
            h_tau(j) = exp(-lambda*tau)*lambda*tau;
            sumh = h_tau(j) +  exp(-lambda*tau);
            
            while sumh < 1 - epsilon
                h_tau = [h_tau; h_tau(j)*lambda*tau/(j+1)];
                j = j + 1;
                sumh = sumh + h_tau(j);
            end
            Pterms = j;
            
            %Compute terms K 
            Kold = cell(1,Pterms);
            Kold{1} = eye(m);
            for j = 2:Pterms
                Kold{j} = Kold{j-1}*P0;
            end
            
            %P(0,tau)
            P0tau = expm(D0*tau);
            Ptot = P0tau;
            
            %P(k,tau), k>0
            k = 1;
            probCum = min(sum(Ptot,2));
            while probCum < 1 - epsilon
                K{1} = P1*Kold{1};
                if k > Pterms
                    htemp = poisspdf(k, lambda*tau);
                else
                    htemp = h_tau(k);
                end
                Ptau{k} = htemp*K{1};
                
                for j = 2:Pterms
                    K{j} = P0*K{j-1} + P1*Kold{j};
                    htemp = htemp*lambda*tau/(k-1+j);
                    Ptau{k} = Ptau{k} + htemp*K{j};
                    j=j+1;
                end
                
                Ptot = Ptot + Ptau{k};

                for j = 1:Pterms
                    Kold{j} = K{j};
                end
                probCum = min(sum(Ptot,2));
                k = k + 1;
            end
            %Ptot
            probCum=sum(Ptot,2);

            %Determine blocks of matrix Ptilde, depending in size of pi
            Kmax = size(Ptau,2);
            numBlocks = size(pi,2)/m;

            %Compute vectors g (h in Neuts)
            sumPi = pi(1:m);
            
            transP0tau = P0tau';
            
            g = zeros(numBlocks,m);
            g(1,:) = linsolve(transP0tau,sumPi')';
            
            i = 2;
            while sum(g(i-1,:),2) < WTacum && i <=numBlocks
                sumPi = sumPi + pi((i-1)*m+1:i*m);
                Kstart = max(0,i-1-Kmax);
                sumG = g(Kstart+1,:)*Ptau{i-1-Kstart};
                for k = Kstart+2:i-1
                    sumG = sumG + g(k,:)*Ptau{i-k};
                end
                g(i,:) = linsolve(transP0tau, (sumPi - sumG)')';
                i=i+1;
            end
            g(i:numBlocks,:)=ones(numBlocks-i+1,1)*thetaA;
            %plot(g)
            
            %Compute Waiting Time distribution at t = {tau,s+tau,2s+tau,...}
            wt = [zeros(1,size(w,2)-1) 1];
            wt(1) = sum(g(c,:));
            for i = 2:size(w,2)-1
                wt(i) = sum(g(i*c,:));
            end
            w = [w;wt];
        end
        w = reshape(w,1,[]);
    end
end