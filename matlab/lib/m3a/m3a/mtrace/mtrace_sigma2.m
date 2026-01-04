function sigma = mtrace_sigma2(T,L)
% Computes the empirical probability of observing a specific 3-element
% sequence of events, i.e. the two-step class transition probabilities.
% Input:
%   T: the inter-arrival times
%   L: the event labels
% Output:
%   sigma: the matrix whose (i,j) element is the probability of observing
%   an event of class i followed by an event of class j.

marks = unique(L);
C = length(marks);

sigma = zeros(C,C,C);

for i = 1:C
    for j = 1:C
        for h = 1:C
            sigma(i,j,h) = sum(L(1:end-2)==marks(i) & ...
                               L(2:end-1)==marks(j) & ...
                               L(3:end)==marks(h));
            sigma(i,j,h) = sigma(i,j,h) / (length(L)-2);
        end
    end
end