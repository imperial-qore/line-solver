%{
%{
 % @file pfqn_joint.m
 % @brief Compute joint queue-length probability distribution.
%}
%}

%{
%{
 % @brief Compute joint queue-length probability distribution.
 % @fn pfqn_joint(n, L, N, Z, lGN)
 % @param n Queue-length state vector (total or per-class).
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector (optional).
 % @param lGN Log normalizing constant at N (optional, computed if not provided).
 % @return pjoint Joint probability of state n.
%}
%}
function pjoint = pfqn_joint(n,L,N,Z,lGN)
% pjoint = pfqn_joint(n,L,N,Z,lGN)
%
% Compute the joint queue-length probability for vector n=(n_1,...,n_M) or
% (n_{11},...n_{M,R}), with M the number of queues, and R the number of
% classes. If there is a think time Z, then n is just the queue population.
% Optionally, the logarithm of the normalizing constant G at N can be passed
% as an input to reduce the computational cost.
%
% Examples:
% * Total queue-lengths
%   pjoint=[];
%   for n=0:4 % queue 1
%       for z=0:4-n % think time
%           pjoint(end+1) = pfqn_joint([n;4-n-z],[10,2;5,4],[2,2],[91,92]);
%       end
%   end
%   sum(pjoint)
%
% * Per-class queue-lengths
%   pjoint=[];
%   for n1=0:4 % queue 1, class 1
%       for n2=0:3  % queue 1, class 2
%           for z1=0:4-n1 % think time, class 1
%               for z2=0:3-n2  % think time, class 2
%                   pjoint(end+1) = pfqn_joint([n1,n2;4-n1-z1,3-n2-z2],[10,2;5,4],[4,3],[91,92]);
%               end
%           end
%       end
%   end
%   sum(pjoint)

[M,R] = size(L);
if nargin<4
    Z= zeros(1,R);
end
if nargin<5
    [~,lGn] = pfqn_ca(L,N,Z);
end
if size(n,2)==1
    % Joint probability of total queue lengths
    if sum(Z)>0
        n0 = sum(N) - sum(n);
        Fjoint = Fper([L;Z],N,[n(:);n0]);
        pjoint = exp(log(Fjoint)-lGn-factln(sum(n0)));
    else
        Fjoint = Fper(L,N,[n(:)]);
        pjoint = exp(log(Fjoint)-lGn);
    end
elseif size(n,2)==R
    n0 = N - sum(n,1);
    % joint probability of total per-class queue-lengths
    if sum(Z)>0
        Fjoint = sum(n0.*log(Z))-sum(factln(n0));
    else
        Fjoint = 0;
    end
    for i=1:M
        Fjoint = Fjoint + multinomialln(n(i,:)) + sum(n(i,:).*log(L(i,:)));
    end
    pjoint = exp(Fjoint-lGn);
else
    line_error(mfilename,'Invalid argument to pfqn_joint');
end
end

function [F,A] = Fper(L,N,m)
[M,R]=size(L);
Ak = [];
for r=1:R
    Ak = [Ak, repmat(L(:,r),1,N(r))];
end
A = [];
for i=1:M
    if m(i)>0
        A((end+1):(end+m(i)),:) = repmat(Ak(i,:),m(i),1);
    end
end
F=perm(A)/prod(factorial(N));
end
