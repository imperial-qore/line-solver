function JM = mtrace_joint(T, A, i)
% Given a multi-class trace, computes the empirical class-dependent joint
% moments that estimate E[ ( X^(a)_j )^i(1) (X^(a)_(j+l) )^i(2) ]
% for all classes a.
%
% Input:
% T: the inter-event times
% A: the class of each event
% i: a vector specifying the power of each variable in the joint moment
%
% Output:
% JM: the joint moment of each class

% number of classes
C = max(A);

% number of events
N = length(A);

% result 
JM = zeros(C,1);

for a = 1:C
    % events of class a, excluding the first and last event
    Ta = T(A(2:(N-1))==a);
    % number of events of class a, excluding the first and last event
    Na = length(Ta);
    tmp = 0;
    for j = 1:(N-2)
       if A(j+1) ==  a
          tmp = tmp + T(j)^i(1) * T(j+1)^i(2);
       end
    end
    JM(a) = 1/Na * tmp;
end

end