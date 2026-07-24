function [AoI_g,AoI_A,AoI_h,AoI_mean,AoI_var,PAoI_g,PAoI_A,PAoI_h,PAoI_mean,PAoI_var]=solveBufferless(tau,T,sigma,S,p)
% solveBufferless calculates parameters for the distributions of AoI and
% pAoI and their first and second moments with the means of MFQ for bufferless systems.
%
% Original source: aoi-fluid toolbox
% Copyright (c) 2020, Ozancan Dogan, Nail Akar, Eray Unsal Atay
% BSD 2-Clause License
%
% Integrated into LINE solver for Age of Information analysis.

% Inputs:
% tau,T   : infitesimal generator and initial probability vector
%           of the arrival process
% sigma,S : infitesimal generator and initial probability vector
%           of the service process
% p       : packet preemption probability

% Outputs:
% AoI_g,AoI_A,AoI_h    : matrix exponential parameters of AoI
% PAoI_g,PAoI_A,PAoI_h : matrix exponential parameters of pAoI
% AoI_mean,AoI_var     : mean and variance of AoI
% PAoI_mean,PAoI_var   : mean and variance of pAoI

k     = size(T,2); % size of the arrival matrix
l     = size(S,2); % size of the service matrix
kappa = -T*ones(k,1);
nu    = -S*ones(l,1);
z     = 2*k*l+k+1; % z is the system size
a     = 1;         % a is the number of negative drift states
b     = z-1;       % b is the number of positive drift states
% construction of the matrix Q in Eqn (13)
Q11 = kron(eye(k),S) + kron(T, eye(l)) + (1-p)*kron(kron(kappa,tau),eye(l));
Q33 = Q11 + p*kron(kappa,kron(ones(l,1),kron(tau,sigma)));
Q   = [Q11,kron(eye(k), nu),zeros(k*l,k*l),p*kron(kappa,ones(l,1)); ...
    zeros(k,k*l),T,kron(kappa,kron(tau,sigma)),zeros(k,1);...
    zeros(k*l,k*l),zeros(k*l,k),Q33,kron(ones(k,1),nu);
    zeros(1,z)];
% construction of the drift matrix R in Eqn (13)
R      = eye(z);
R(z,z) = -1;
% Qtilde Construction in Eqn (13)
Qtilde          = Q;
Qtilde(z,1:k*l) = kron(tau,sigma);
Qtilde(z,z)     = -1;
% construction of the orthogonal matrix P in Line 2 of Alg. 1
QR = Q/R;
u1 = [ones(z-1,1);-1];
u2 = [1;zeros(z-1,1)];
u  = u1-norm(u1,2)*u2;
P  = eye(z)-(2*(u)*(u'))/(u'*u);
% construction of the matrices A, H, and Qtildestar in Line 3 of Alg. 1
Ae         = P'*QR*P;
A          = Ae(a+1:z,a+1:z);
H          = P(1:z,a+1:z)';
Qtildestar = Qtilde(b+1:z,1:z);
%solve for the vectors g and d in Line 4 of Alg. 1
EqnMatrix     = [H*R,-inv(A)*H*ones(z,1);-Qtildestar ones(a,1)];
RightHandSide = [zeros(1,z) 1];
Solution      = mldivide(EqnMatrix',RightHandSide')';
g             = Solution(1:b);
d             = Solution(b+1:z);
% We now have the ME representation in Eqn (8) for the MFQ X(t)
% For AoI, restrict to states in Phases 2 and 3 only
SelectedStates = [zeros(k*l,1);ones(k*l+k,1);0];
AoI_h          = H*SelectedStates;
NormConstant   = -g/A*AoI_h;
AoI_g          = g/NormConstant;
AoI_A          = A;
% now we have the general form in Eqn (4)for the AoI distribution
AoI_mean = AoI_g/(AoI_A^2)*AoI_h; % the moment expressions in Eqn (4)
AoI_var  = -2*AoI_g/(AoI_A^3)*AoI_h-AoI_mean^2;
% For PAoI, restrict to states in Phase 3 scaled suitably by the entries of nu
SelectedStates = [zeros(k*l+k,1);kron(ones(k,1),nu);0];
PAoI_h         = H*SelectedStates;
NormConstant   = -g/A*PAoI_h;
PAoI_g         = g/NormConstant;
PAoI_A         = A;
% now we have the general form in Eqn (4)for the AoI distribution
PAoI_mean = PAoI_g/(PAoI_A^2)*PAoI_h; % the moment expressions again in Eqn (4)
PAoI_var  = -2*PAoI_g/(PAoI_A^3)*PAoI_h-PAoI_mean^2;
end
