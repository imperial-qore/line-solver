function demandEst = des_FMLPS(rt,class,ql,lambda,V,W)
% DES_FMLPS implements the FMLPS demand estimation method
% RT:       response times samples. 
%           column vector with all response time samples  
% CLASS:    class of the request samples
%           column vector with the class of each sample
% QL:       queue length samples
%           matrix with R columns containing the number of jobs 
%           of each class observed by each sample
% LAMBDA:   estimated think time rate
%           1xR vector with the think time for each class as entries
% V:        number of cores
% W:        number of threads  
%
% Copyright (c) 2012-2014, Imperial College London 
% All rights reserved.


R = length(lambda); % number of customer types

xLB = min(rt)*ones(1,R)/W;        % lower bounds on x variables
xUB = max(rt)*ones(1,R);  % upper bounds on x variables

% initial point
meanQL = mean(sum(ql,2));
Vtilde = min(meanQL,V);
x0 = zeros(1,R);
for j = 1:R
    if ~isempty(rt(class==j))
        x0(j) = Vtilde * mean(rt(class==j)) / meanQL;
    else
        x0(j) = xLB(j);
    end
end

%% options
options = optimset();
options.Display = 'iter';
options.LargeScale = 'on';
options.MaxIter =  1e10;
options.MaxFunEvals = 1e10;
options.MaxSQPIter = 5000;
options.TolCon = 1e-6;
options.Algorithm = 'interior-point';
options.TolX = 1e-10;


[demandEst, fObjFun]=fmincon(@objfun,x0,[],[],[],[],xLB,xUB,[],options);
fObjFun = - fObjFun; % change sign

%% objective (likelihood) function
    function f = objfun(x)        
        TOL = 1e-6;
        rates = 1./x; % x = mean demands
        ftemp = zeros(size(rt));
        
        %create CQN object
        sched = {'inf'; 'ps'};
        K = size(ql,2);
        P = zeros(2*K, 2*K);
        P(1:K,K+1:2*K) = eye(K);
        P(K+1:2*K, 1:K) = eye(K);
        fullRates = [lambda; rates];
        myCQN = CMCQNCS(2, K, 0, [-1; V], fullRates, sched,P, zeros(1,K));
        
        
        parfor r = 1:length(rt)
        %for r = 1:length(rt)
            ftemp(r) = log(TOL + des_FMLPS_like(myCQN, rt(r), class(r), ql(r,:), W) );
        end
        f = -sum(ftemp);
    end

end