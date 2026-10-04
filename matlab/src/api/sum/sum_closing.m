function [XN,QN,UN,RN,TN,it] = sum_closing(lambda0,scva,L,mi,scv,N,Z,Kclosed,tol,maxiter)
% [XN,QN,UN,RN,TN,IT] = SUM_CLOSING(LAMBDA0,SCVA,L,MI,SCV,N,Z,KCLOSED,TOL,MAXITER)
%
% Closing method for open and mixed non-product-form queueing networks
% (Bolch et al., Sec. 10.1.5), solved with the summation method.
%
% The external world of each open class is replaced by an additional
% -/G/1 station with service rate mu_inf,r = Ropen*lambda0(r), where Ropen
% is the number of open classes, service SCV equal to the interarrival
% time SCV of the open class, and unit visit ratio. The resulting closed
% network is then solved by SUM_CLOSED with a large population KCLOSED
% for the open classes (default: 5000, as recommended for the summation
% method). Closed classes are passed through unchanged, which makes the
% method applicable to mixed networks.
%
% Input:
% lambda0 - 1xR external arrival rates (0 for closed classes)
% scva    - 1xR interarrival time SCVs of the open classes (1 if Poisson)
% L       - MxR service demand matrix of the original network, with visit
%           ratios of open classes normalized per external arrival
% mi      - Mx1 number of servers (Inf for infinite-server stations)
% scv     - MxR service time SCVs (pass 1 for insensitive stations)
% N       - 1xR populations: Inf (or NaN) for open classes, finite
%           integers for closed classes
% Z       - 1xR think times
% Kclosed - closing population for the open classes (default: 5000)
% tol     - convergence tolerance (default: 1e-6)
% maxiter - maximum number of iterations (default: 10000)
%
% Output:
% XN - 1xR class throughputs (for open classes, XN approaches lambda0
%      from below as KCLOSED grows)
% QN - MxR mean queue lengths at the original stations
% UN - MxR utilizations at the original stations
% RN - MxR residence times at the original stations
% TN - 1xR mean response time in the original network, TN=sum(QN)./XN
% it - number of iterations
%
% Reference: G. Bolch, S. Greiner, H. de Meer, K.S. Trivedi, Queueing
% Networks and Markov Chains, 2nd ed., Wiley, 2006, Sec. 10.1.5.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

[M,R] = size(L);
lambda0 = lambda0(:)';
if nargin<2 || isempty(scva), scva = ones(1,R); end
if nargin<4 || isempty(mi), mi = ones(M,1); end
if nargin<5 || isempty(scv), scv = ones(M,R); end
if nargin<6 || isempty(N), N = Inf*ones(1,R); end
if nargin<7 || isempty(Z), Z = zeros(1,R); end
if nargin<8 || isempty(Kclosed), Kclosed = 5000; end
if nargin<9 || isempty(tol), tol = 1e-6; end
if nargin<10 || isempty(maxiter), maxiter = 10000; end
scva = scva(:)';
N = N(:)';
Z = Z(:)';

openClasses = find(lambda0>0);
Ropen = length(openClasses);
if Ropen==0
    line_error(mfilename,'No open class: use sum_closed for closed networks.');
end

% augment with the closing -/G/1 station, visited by open classes only
Laug = [L; zeros(1,R)];
scvaug = [scv; ones(1,R)];
Naug = N;
for r=openClasses
    Laug(M+1,r) = 1/(Ropen*lambda0(r));
    scvaug(M+1,r) = scva(r);
    Naug(r) = Kclosed;
end
miaug = [mi(:); 1];

[XN,QNa,UNa,RNa,it] = sum_closed(Laug,Naug,Z,miaug,scvaug,tol,maxiter);

QN = QNa(1:M,:);
UN = UNa(1:M,:);
RN = RNa(1:M,:);
TN = zeros(1,R);
for r=1:R
    if XN(r)>0
        TN(r) = sum(QN(:,r))/XN(r);
    end
end
end
