function [ pi0, En1 ] = computePi(T, arrival, services, service_h, C, S, A_jump)
%
% computePi(T, arrival, service, C) returns the steady state distribution
% for the case with 2 replicas (2 seperate queues)
%

if services.SerChoice  == 1
    
    %% Construct the Se, Sestar matrices
    ms = size(S,2);
    [S_notallbusy] = constructNotAllBusy(C, services, service_h);
    Q0 = kronsum(S_notallbusy,arrival.lambda0);
    
    da = size(arrival.lambda0,1);
    Igral = lyap(T,kron(eye(ms),arrival.lambda0),-eye(da*ms));
    pi0mat = Igral*kron(A_jump,eye(da))*inv(Q0)*kron(eye(ms),arrival.lambda1);
    
    % compute pi0 = pi0*pi0mat
    pi0 = [zeros(1,ms*da),-1]/[(pi0mat - eye(ms*da)),sum(inv(T),2)];
    
    En1 = 1/sum(pi0)*sum(pi0*pi0mat); % Igral*kron(A_jump,eye(da))*inv(Q0)*kron(eye(ms),arrival.lambda1))
    
else
    
    [Se, Sestar, R0, Ke, Kc] = constructSRK(C, services, service_h, S);
    dtmat = size(T,1)/size(arrival.lambda0,2);
    dsexp = size(Se,1);
    
    Sedash = Se((dtmat+1):dsexp,(dtmat+1):dsexp);
    Rbusy = R0((dtmat+1):dsexp,:);
    Iidle = eye(dsexp);
    Iidle = Iidle(:,(dtmat+1):dsexp);
    da = length(arrival.lambda0);
    Iidle = kron(Iidle,eye(da));
    
    Q0 = kronsum(Sedash,arrival.lambda0); % state changes between arrivals
    Qbusy = kron(Rbusy,arrival.lambda1); % return to all-busy state
    Qidle = Q0;
    Bmap = -(Iidle/Qidle)*Qbusy;
    Kemap = kron(Ke,eye(da));
    Kcmap = kron(Kc,eye(da));
    
    BB = kron(eye(size(Sestar,2)),arrival.lambda0);
    Igral = lyap(T,BB,Kemap*kron(Sestar,eye(da))); % solves the Sylvester equation AX+XB+C=0
    pi0mat = Igral*Bmap*Kcmap; % equation (7.2)
    
    % compute pi0 = pi0*pi0mat
    pi0 = [zeros(1,dtmat*da),-1]/[(pi0mat - eye(dtmat*da)),sum(inv(T),2)];
    
    En1 = 1/sum(pi0)*sum(-pi0*Igral*Iidle/Qidle*kron(eye(size(Sedash,2)),arrival.lambda1)); % E(n1)
    
end