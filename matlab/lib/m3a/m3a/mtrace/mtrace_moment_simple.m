function MC = mtrace_moment_simple(T,L,k)
% Computes the k-th order moment of the inter-arrival time between an event
% of class i and an event of class j, for all possible pairs of classes.
% Input
%   T: inter-arrival times
%   L: labels
%   k: order of the moment
% Output
%   MC: the element in (i,j) is the k-th order moment of the inter-arrival
%       time between an event of class i and an event of class j

marks = unique(L);
C = length(marks);

MC = zeros(C,C);
count = zeros(C,C);

for t = 2:length(T)
    for i = 1:C
        for j = 1:C
            if L(t-1) == marks(i) && L(t) == marks(j)
                MC(i,j) = MC(i,j) + T(t)^k;
                count(i,j) = count(i,j) + 1;
            end
        end
    end
end

MC = MC ./ count;