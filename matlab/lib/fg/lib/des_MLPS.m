function demandEst = des_MLPS(rt,class,ql,lambda,nCores)
% DES_MLPS implements the MLPS demand estimation method
% RT:       response times samples. 
%           column vector with all response time samples  
% CLASS:    class of the request samples
%           column vector with the class of each sample
% QL:       queue length samples
%           matrix with R columns containing the number of jobs 
%           of each class observed by each sample
% LAMBDA:   estimated think time rate
%           1xR vector with the think time for each class as entries
% NCORES:   number of cores
%
% Copyright (c) 2012-2014, Imperial College London 
% All rights reserved.

R = length(lambda); % number of customer types
% initial point
meanQL = mean(sum(ql,2));
Vtilde = min(meanQL,nCores);
x0 = zeros(1,R);
for j = 1:R
    x0(j) = Vtilde * mean(rt(class==j)) / meanQL;
end

xLB = zeros(1,R);        % lower bounds on x variables
xUB = max(rt)*ones(1,R);  % upper bounds on x variables


%% options
options = optimset();
options.Display = 'iter';
options.LargeScale = 'on';
options.MaxIter =  1e10;
options.MaxFunEvals = 1e10;
options.MaxSQPIter = 5000;
options.TolCon = 1e-6;
options.Algorithm = 'interior-point';

[demandEst, fObjFun]=fmincon(@objfun,x0,[],[],[],[],xLB,xUB,[],options);
fObjFun = - fObjFun; % change sign


%% objective (likelihood) function
    function f = objfun(x)        
        TOL = 1e-6;
        rates = 1./x; % x = mean demands
        ftemp = zeros(size(rt));
        parfor r = 1:length(rt)
            ftemp(r) = log(TOL + des_MLPS_like(rt(r),class(r),ql(r,:),rates,lambda,nCores));
        end
        f = -sum(ftemp);
    end

end