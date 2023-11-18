function LIKE = des_FMLPS_like(myCQN, Rsam, tagClass, aQueue, W)
% DES_FMLPS_LIKE computes the likelihood of the response time Rsam based on
% a fluid model
% myCQN:    queueing network model
% Rsam:     observed response time for tagged job 
% tagClass: class of tagged job
% aQueue:   number of jobs of each class seen in the system upon arrival (including arriving job)
%           1xR vector of integers            
% W:        number of threads/jobs in the model
%
% Copyright (c) 2012-2014, Imperial College London 
% All rights reserved.

% estimate number of jobs in the delay station 
delayJobs = (W-sum(aQueue))*myCQN.rates(1,:)/sum(myCQN.rates(1,:)); 
NK = aQueue + delayJobs; 
N = W;

myCQN.NK = NK;
myCQN.N = N;

myCQNCSCox = CMCQNCS2CMCQNCSCox(myCQN);
y0 = [delayJobs aQueue];
LIKE = CMCQN_CS_fluid_ps_RT_like(myCQNCSCox, y0, tagClass, Rsam);


