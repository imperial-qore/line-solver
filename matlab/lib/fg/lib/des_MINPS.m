function demandEst = des_MINPS(rt,class,ql,lambda,V)
% DES_MINPS implements the MINPS demand estimation method
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
%
% Copyright (c) 2012-2014, Imperial College London 
% All rights reserved.

% number of classes
R = max(class); 

% put samples in RPS format
rtRPS = cell(1,R);
qlRPS = cell(1,R);
for r = 1:R
    rtRPS{r} = rt(class==r); 
    qlRPS{r} = ql(class==r,:); 
end


% MINPS
% run MLPS
demandEstMLPS = des_MLPS(rt,class,ql,lambda,V);
% run RPS
demandEstRPS = des_RPS(rt,class,ql,V);
% choose the smallest mean result
if mean(demandEstMLPS) < mean(demandEstRPS)
    demandEst = demandEstMLPS;
else
    demandEst = demandEstRPS;
end