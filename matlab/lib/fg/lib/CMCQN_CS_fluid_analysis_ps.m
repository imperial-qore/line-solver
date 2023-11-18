function [ode_h,ode_rates_h,ode_jumps_h,q_indices] = CMCQN_CS_fluid_analyzer_ps(N, Mu, Phi, P, S)
% CMCQN_CS_FLUID_ANALYZER_PS provides the methods to perform the ODE 
% analysis of a Closed Multi-Class Queueing Network with Class Switching
% with Coxian Service Times (CMCQNCSCox). 
% More details on this type of queueing networks can be found 
% on the LINE documentation, available at http://code.google.com/p/line
%
% Parameters:
% N:    total number of jobs
% Mu:   service rates in each station (in each phase), for each stage
% Phi:  probability of service completion in each stage of the service 
%       process in each station, for each stage
% P:    routing matrix for each stage
% S:    number of servers in each station, for each stage
%
% Output:
% ode_h:        handler of the ODE system
% q_indices:    indices of each job class and station, in the state vector
% ode_jumps_h:  handler of the jumps in the ODE system
% ode_rates_h:  handler of the transition rates in the ODE system
%
% Copyright (c) 2012-2014, Imperial College London 
% All rights reserved.
        
    
M = length(S);    % number of stations
K = size(Mu,2);   % number of classes

match = zeros(M,K); % indicates whether a class is served at a station
for i = 1:M
    for c = 1:K
        match(i,c) = sum( P(:,(i-1)*K+c) ) > 0 || sum( P((i-1)*K+c,:) )>0;
    end
end

function xj = get_index(j,k)
    % n is the state vector
    % j is the queue station index
    % k is the class index
    % returns the index of the queue-length element xi! in the state description
    xj = 1;
    for z = 1 : (j-1)
        for y = 1:K
            xj = xj + Kic(z,y);
        end
    end
    for y = 1:k-1
        xj = xj + Kic(j,y);
    end

end

function rates = ode_rates(x)
    rates = [];
    n = zeros(1,M); % total number of jobs in each station
    for i = 1:M
        for c = 1:K
            xic = q_indices(i,c);
            n(i) = n(i) + sum( x(xic:xic+Kic(i,c)-1 ) );
        end
        if S(i) == sum(N)
            n(i) = 1;
        end
    end



    for i = 1 : M           %transition rates for departures from any station to any other station
        for c = 1:K         %considers only transitions from the first service phase (enough for exp servers)
            if match(i,c)>0
            xic = q_indices(i,c);
            for j = 1 : M
                for l = 1:K
                    if P((i-1)*K+c,(j-1)*K+l) > 0 
                    for k = 1:Kic(i,c)
                        %pure ps + fcfs correction
                        if x(xic+k-1) > 0 && n(i) > S(i)
                            rates = [rates; Phi{i,c}(k) * P((i-1)*K+c,(j-1)*K+l) * Mu{i,c}(k) * x(xic+k-1)/n(i) * S(i);]; % f_k^{dep}
                        elseif x(xic+k-1) > 0
                            rates = [rates; Phi{i,c}(k) * P((i-1)*K+c,(j-1)*K+l) * Mu{i,c}(k) * x(xic+k-1);]; % f_k^{dep}
                        else
                            rates = [rates; 0;]; % f_k^{dep}
                        end
                    end
                    end
                end
            end
            end
        end
    end

    for i = 1 : M           %transition rates for "next service phase" (phases 2...)
        for c = 1:K
            if match(i,c)>0
            xic = q_indices(i,c);
            for k = 1 : (Kic(i,c) - 1)
                if x(xic+k-1) > 0 
                    rates = [rates; (1-Phi{i,c}(k))*Mu{i,c}(k)*x(xic+k-1)/n(i)];
                else
                    rates = [rates; 0 ]
                end
            end
            end
        end
    end

end

function d = ode_jumps(x)
    %reshape(x,K,M)'
    d = [];         %returns state changes triggered by all the events
    for i = 1 : M   %state changes from departures in service phases 2...
        for c = 1:K
            if match(i,c)>0
            xic = q_indices(i,c);
            for j = 1 : M
                for l = 1:K
                    if P((i-1)*K+c,(j-1)*K+l) > 0 
                    xjl = q_indices(j,l);
                    for k = 1 : Kic(i,c)
                        jump = zeros(length(x),1);
                        jump(xic) = jump(xic) - 1; %type c in stat i completes service
                        jump(xjl) = jump(xjl) + 1; %type c job starts in stat j
                        d = [d jump;];
                    end
                    end
                end
            end
            end
        end
    end

    for i = 1 : M   %state changes from "next service phase" transition in phases 2...
        for c = 1:K
            if match(i,c)>0
            xic = q_indices(i,c);
            for k = 1 : (Kic(i,c) - 1)
                jump = zeros(length(x),1);
                jump(xic+k-1) = jump(xic+k-1) - 1;
                jump(xic+k) = jump(xic+k) + 1;
                d = [d jump;];
            end
            end
        end
    end        
end


function diff = ode(t,x)
    diff = ode_jumps(x)*ode_rates(x);%rate of change in state x
end

ode_h = @ode;
ode_rates_h = @ode_rates;
ode_jumps_h = @ode_jumps;

q_indices = zeros(M,K); 
Kic = zeros(M,K);
for i = 1 : M
    for c = 1:K
        Kic(i,c) = length(Mu{i,c});
        q_indices(i,c) = get_index(i,c);
    end
end

end
