function sigma = mtrace_sigma(T,L)
% Computes the empirical probability of observing a specific 2-element
% sequence of events, i.e. the one-step class transition probabilities.
% Input:
%   T: the inter-arrival times
%   L: the event labels
% Output:
%   sigma: the matrix whose (i,j) element is the probability of observing
%   an event of class i followed by an event of class j.

marks = unique(L);
C = length(marks);

sigma = zeros(C,C);

for i = 1:C
    for j = 1:C
        sigma(i,j) = sum(L(1:end-1)==marks(i) & L(2:end)==marks(j));
        sigma(i,j) = sigma(i,j) / (length(L)-1);
    end
end