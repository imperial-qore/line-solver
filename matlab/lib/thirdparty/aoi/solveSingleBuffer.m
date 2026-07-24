function [AoI_g,AoI_A,AoI_h,AoI_mean,AoI_var,PAoI_g,PAoI_A,PAoI_h,PAoI_mean,PAoI_var] = solveSingleBuffer(lambda,sigma,S,r)
% solveSingleBuffer calculates parameters for the distributions of AoI and
% pAoI and their first and second moments with the means of MFQ for single buffer systems.
%
% Original source: aoi-fluid toolbox
% Copyright (c) 2020, Ozancan Dogan, Nail Akar, Eray Unsal Atay
% BSD 2-Clause License
%
% Integrated into LINE solver for Age of Information analysis.

% Inputs:
% lambda  : arrival rate
% sigma,S : infitesimal generator and initial probability vector
%           of the service process
% r       : packet replacement probability
% Outputs:
% AoI_g,AoI_A,AoI_h    : matrix exponential parameters of AoI
% PAoI_g,PAoI_A,PAoI_h : matrix exponential parameters of pAoI
% AoI_mean,AoI_var     : mean and variance of AoI
% PAoI_mean,PAoI_var   : mean and variance of pAoI

l  = size(S,2);    % size of the service matrix
nu = -S*ones(l,1);
z  = l+2;          % z is the system size of waiting matrix
a  = 2;            % a is the number of negative drift states
b  = l;            % b is the number of positive drift states
% construction of the matrix Q in Eqn (16)
Q = [S,nu,zeros(l,1);
    zeros(1,l),-lambda,lambda;
    zeros(1,l+2)];
% construction of the drift matrix R in Eqn (16)
R          = eye(z);
R(z-1,z-1) = -1;
R(z,z)     = -1;
% Qtilde Construction in Eqn (16)
Qtilde=[zeros(l,l+2);
    lambda*sigma,-lambda,0;
    sigma,0,-1];
% Steady state distribution of MFQ
e   = ones(z,1);
pik = e'/(Q+e*e');
% construction of the orthogonal matrix P
% in Line 2 of Alg. 1 with Schur Algorithm
QR     = Q/R;
x_R    = R*e;
x_L    = pik;
A1     = QR+(x_R*x_L)/(x_L*x_R);
[U,TT] = schur(A1);
[P,~]  = ordschur(U,TT,'rhp');
% construction of the matrices A, H, and Qtildestar in Line 3 of Alg. 1
Ae         = P'*QR*P;
A          = Ae(a+1:z,a+1:z);
H          = P(1:z,a+1:z)';
Qtildestar = Qtilde(b+1:z,1:z);
%solve for the vectors g and d in Line 4 of Alg. 1
EqnMatrix     = [H*R,-inv(A)*H*ones(z,1);-Qtildestar ones(a,1)];
RightHandSide = [zeros(1,z) 1];
Solution      = mldivide(EqnMatrix',RightHandSide')';
wait_g        = Solution(1:b);
wait_d        = Solution(b+1:z);
% We now have the ME representation in Eqn (8) for the MFQ X(t)
c_0    = wait_d(1);
wait_A = A-r*lambda*eye(l);
wait_H = H*[zeros(l,1);1;r];
n_1    = 1/(-wait_g/wait_A*wait_H + c_0);
wait_g = n_1*wait_g;
% Lemma 1
Mdiag  = -wait_A\wait_H;
M      = diag(Mdiag);
B      = M\wait_A*M;
beta   = wait_g*M;
beta_0 = 1-sum(beta);
psi    = -B*ones(l,1);
% construction of the matrix Q in Eqn (19)
z = 4*l+2; % z is the system size
a = 1;     % a is the number of negative drift states
b = z-1;   % b is the number of positive drift states

% construction of the matrix Q in Eqn (19)
Q = [B,kron(psi,sigma),zeros(l,2*l+2);
    zeros(l,l),S-lambda*eye(l),lambda*eye(l),nu,zeros(l,l+1);
    zeros(l,2*l),S,zeros(l,1),kron(nu,sigma),zeros(l,1);
    zeros(1,3*l),-lambda,lambda*sigma,zeros(1,1);
    zeros(l,3*l+1),S,nu;
    zeros(1,4*l+2)];
% construction of the drift matrix R in Eqn (19)
R      = eye(z);
R(z,z) = -1;
% Qtilde Construction in Eqn (19)
Qtilde            = Q;
Qtilde(z,1:l)     = beta;
Qtilde(z,l+1:2*l) = beta_0*sigma;
Qtilde(z,z)       = -1;
% construction of the orthogonal matrix P in Line 2 of Alg. 1
QR = Q/R;
u1 = [ones(z-1,1);-1];
u2 = [1;zeros(z-1,1)];
u  = u1-norm(u1,2)*u2;
P  = eye(z)-(2*(u)*(u'))/(u'*u);
% construction of the matrices A, H, and Qtildestar in Line 3 of Alg. 1
Ae         = P'*(QR)*P;
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
% For AoI, restrict to states in Phases 4 and 5 only
SelectedStates = [zeros(3*l,1);ones(l+1,1);0];
AoI_h          = H*SelectedStates;
NormConstant   = -g/A*AoI_h;
AoI_g          = g/NormConstant;
AoI_A          = A;
% now we have the general form in Eqn (4)for the AoI distribution
AoI_mean = AoI_g/(AoI_A^2)*AoI_h; % the moment expressions in Eqn (4)
AoI_var  = -2*AoI_g/(AoI_A^3)*AoI_h-AoI_mean^2;
% For PAoI, restrict to states in Phase 5 scaled suitably by the entries of nu
SelectedStates = [zeros(3*l+1,1);nu;0];
PAoI_h         = H*SelectedStates;
NormConstant   = -g/A*PAoI_h;
PAoI_g         = g/NormConstant;
PAoI_A         = A;
% now we have the general form in Eqn (4)for the AoI distribution
PAoI_mean = PAoI_g/(PAoI_A^2)*PAoI_h; % the moment expressions again in Eqn (4)
PAoI_var  = -2*PAoI_g/(PAoI_A^3)*PAoI_h-PAoI_mean^2;
end
