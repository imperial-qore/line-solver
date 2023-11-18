function LIKE = des_MLPS_like(Rsam,tagClass, aQueue, rates, muZ, nCores)
% DES_MLPS_LIKE 
% Rsam:     observed response time for tagged job 
% tagClass: class of tagged job
% aQueue:   number of jobs of each class seen in the system upon arrival (including arriving job)
%           1xR vector of integers            
% rates:    exponential service rates of each class 
%           1xR vector of doubles
% muZ:      think-time rates
%           1xR vector of doubles    
% nCores:   number of cpu cores
%
% Copyright (c) 2012-2014, Imperial College London 
% All rights reserved.


M = 2;   %number of nodes in the net (one processing, one thinking)
R = length(rates);
% An additional customer class is created to represent the tagged customer.
% The appropriate rate and number are assigned to the vectors N, rates and muZ 
% as well as the routing matrix in the cell P
rates(R+1)=rates(tagClass);
muZ(R+1)=muZ(tagClass);
P = {}; for r=1:R, P{r} = circul(M); end; P{R+1} = P{tagClass};
N = aQueue; N(tagClass)=N(tagClass)-1; N(R+1)=1; % tagged class

f = 1; % station to filter (queue)
[MAPQ,SS] = pfqn_multi_filtration_V2(rates, N, P, muZ, f, R+1,nCores);    

subset = find(SS(:,size(SS,2)-M+1)==1); %considers states where the 
                                        %first station (processing) has
                                        %exactly one job of the tagged
                                        %class - these are the states that
                                        %matter since these are the states
                                        %that the chain visits while the
                                        %tagged job is still in the
                                        %processing queue

A = MAPQ{1};
A = full(A(subset,subset)); %subgenerator matrix for states in subset
SS = SS(subset,1:M:end); %for states in subset, contains representation 
                         %for number of jobs in the first node (processing)
                         %of each type - simplifies representation of the
                         %state space to consider only the processing
                         %station
pie = zeros(1,length(A));
pie(matchrow(SS,N)) = 1; %pie is the initial distribution vector, contains 
                         %a 1 in the position that corresponds to the state
                         %observed upon arrival (in vector N)
MAP = {A,-A*ones(size(pie))'*pie};% builds a MAP with A as generator matrix
                                  % and a renewal (PH)
LIKE = map_pdf(MAP,Rsam);

end