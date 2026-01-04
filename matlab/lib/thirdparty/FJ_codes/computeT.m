function [ T, S, A_jump, S_Arr, sum_Ajump ] = computeT(arrival, services, service_h, C, T_Mode)

% The function calculates T matrix of Equation (2) using both methods in
% Section 5.1: the Sylvester-equation approach (Sylv) and
% 5.2: A Riccati-equation approach (NARE): (default: 'NARE')
%   'Sylves' : solves a Sylvester matrix equation at each step
%                          using a Hessenberg algorithm
%   'NARE  ' : solves the Sylvester matrix equation by defining
%            a non-symmetric algebraic Riccati equation (NARE)
%
%

% arrival is a structure
% that must contain arrival.lambda0 and arrival.lambda1
% where lambda0 is the intensity-matrix for state changes
% in the arrival process not accompanied by arrivals and
% lambda1 is the intensity-matrix for state changes
% accompanied by arrivals.
%
% service is a structure array where the number of structures equals
% the number of servers and a structure contains beta and S
% where service(i).beta and service(i).S corresponds to the representation
% (beta,S) for the servicetime of the i'th server.

% S: without a new Fj job starting service;
% A_jump: with a new FJ job starting service;

[ S, A_jump ] = build_SA( services, service_h, C );

d0 = size(arrival.lambda0,1);
S_Arr = kron(S,eye(d0));
A_jump_Arr = kron(A_jump,speye(d0));

if  (strfind(T_Mode,'Sylvest')>0)
    % The Sylvester-equation approach (Sylv)
    tic
    ms = size(S,1);
    m = ms*d0;
    Tnew = S_Arr; % withour service completion
    Told = zeros(m,m);
    ID0 = kron(eye(ms),arrival.lambda0);
    DS = kron(eye(ms),arrival.lambda1)*A_jump_Arr;
    [U,Tr] = schur(ID0,'complex');
    % U'*Tr*U = A
    while(max(max(abs(Told-Tnew)))>10^(-10))
        Told = Tnew;
        L = Q_Sylvest(U,Tr,Tnew);
        Tnew = S_Arr + L*DS;
    end
    T = Tnew;
    res_norm=norm(T*L+L*kron(eye(ms),arrival.lambda0)+eye(m),inf);
    fprintf('Final Residual Error for T matrix: %d\n',res_norm);
    time_Sylv = toc
else % default
    % NARE
    tic
    T = computeT_NARE(arrival.lambda0, arrival.lambda1, S_Arr, A_jump);
    time_NARE= toc
end

sum_Ajump = sum(A_jump_Arr,2);




