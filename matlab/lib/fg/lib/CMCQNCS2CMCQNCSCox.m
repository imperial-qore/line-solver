function myCQNCox = CMCQNCS2CMCQNCSCox(myCQN)
% B = CMCQNCS2CMCQNCSCOX(A) transforms a CMCQNCS model A (exponential service times) 
% into its equivalent representation as a CMCQNCSCOX model B (Coxian service times) 
%
% Copyright (c) 2012-2014, Imperial College London 
% All rights reserved.


mu = cell(myCQN.M, myCQN.K);
phi = cell(myCQN.M, myCQN.K);
for i =1:myCQN.M
    for k = 1:myCQN.K
        mu{i,k} = myCQN.rates(i,k);
        phi{i,k} = 1;
    end
end
S = myCQN.S;

myCQNCox = CMCQNCSCox(myCQN.M, myCQN.K, myCQN.N, S, mu, phi, myCQN.sched, myCQN.P, myCQN.NK, myCQN.nodeNames,myCQN.classNames);
