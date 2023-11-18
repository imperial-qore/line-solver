function LIKE = CMCQN_CS_fluid_ps_RT_like(myCQN, y0, taggedClass, Rsampled)
% CMCQN_CS_FLUID_PS_RT_LIKE computes the likelihood of the response time
% observed based on the fluid model
% myCQN:        CMCQNCSCox object describing the network
% max_iter:     maximum number of iterations
% y_0:          initial state
% taggedClass:  tagged class
% Rsampled:     response time sampled
%
% Copyright (c) 2012-2014, Imperial College London 
% All rights reserved.



M = myCQN.M;    %number of stations
K = myCQN.K;    %number of classes
N = myCQN.N;    %population
Lambda = myCQN.mu;
Pi = myCQN.phi;
P = myCQN.P;
S = myCQN.S;

refNode = 2;
delayNode = -1;
for i = 1:M
    %Set number of servers in delay station = population
    if S(i) == -1;
        S(i) = N;
        delayNode = i;
    end
end


%% initialization
phases = zeros(M,K);
for i = 1:M;
    for k = 1:K
        phases(i,k) = length(Lambda{i,k});
    end
end

%% response time analysis - starting from y0
% add artifial/tagged class
k = taggedClass;
newK = K + 1;
numNewK = newK - K;
newLambda = cell(M,newK);
newPi = cell(M,newK);
newP = zeros(M*newK, M*newK);

%set lambda and pi for new class
newLambda(:,1:K) = Lambda(:,:);
for l = 1:numNewK
    newLambda(:,K+l) = Lambda(:,k);
end
newPi(:,1:K) = Pi(:,:);
for l = 1:numNewK
    newPi(:,K+l) = Pi(:,k);
end

% extend matrix P to consider tagged class
for l = 1:K
    for m = 1:K
        newP(l:newK:end,m:newK:end) = P(l:K:end,m:K:end);
    end
end
PextraClass = P(taggedClass:K:end,taggedClass:K:end);
PextraClass(refNode,:) = zeros(1,M);
newP(K+1:newK:end,K+1:newK:end) = PextraClass;
newP((refNode-1)*newK+K+1, (delayNode-1)*newK+taggedClass) = 1;

% build ODE model
[newOde_h, newOde_rates_h, newOde_jumps_h] = CMCQN_CS_fluid_analyzer_ps(N, reshape({newLambda{:,:}},M,newK), reshape({newPi{:,:}},M,newK), newP, S);

newPhases = zeros(M,newK);
newPhases(:,1:K) = phases;
for l = K+1:newK
    newPhases(:,l) = phases(:,k);
end

%determine amount of fluid for the new class
newFluid = 1;

%move fluid from original class to the new one
addY = zeros(1, sum(sum(newPhases(:,:))) );
addY( sum(sum(newPhases(1:refNode-1,:))) + k ) = - newFluid; %ref node, class k
addY( sum(sum(newPhases(1:refNode-1,:))) + K+1 ) = newFluid; %ref node, new class
newY0 = zeros(1, sum(sum(newPhases(:,:))));
for i = 1:M
    for l = 1:K
        idx = sum(sum(newPhases(1:i-1,:))) + sum(newPhases(i,1:l-1));
        idxOld = sum(sum(phases(1:i-1,:))) + sum(phases(i,1:l-1));
        newY0( idx + 1: idx + newPhases(i,l)  ) = y0(idxOld+1:idxOld + phases(i,l));
    end
end
newY0 = newY0 + addY;
y0 = newY0;


% solve ODE (fluid model) from 0 to Rsampled
opt = odeset('AbsTol', 1e-8, 'RelTol', 1e-5, 'NonNegative', 1:length(y0),'Events',@events);
[t, yt] = ode15s(newOde_h, [0 Rsampled], y0, opt);

% determine likelihood based on terminal state
if Rsampled <= t(end)
    idx = sum(sum(newPhases(1:refNode-1,: ) )) + newK ; 
    lastState = yt(end,:);
    lastRates = newOde_jumps_h(lastState)*newOde_rates_h(lastState);

    LIKE = -lastRates(idx)/newFluid;
else
    LIKE = 0;
end

return


function [value,isterminal,direction] = events(t,y)
    idx = sum(sum(newPhases(1:refNode-1,:))) + newK;
    value = y(idx);
    isterminal = 1;
    direction = 0;
end


end

