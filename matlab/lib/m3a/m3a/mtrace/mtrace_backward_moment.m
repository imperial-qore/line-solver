function [M] = mtrace_backward_moment(T,A,ORDERS,NORM)
% Computes the backward moments of a marked trace.
% Input:
% - T:      the inter-arrival times
% - C:      the class labels
% - ORDERS: vector with the orders of the moments to compute
% - NORM:   0 to return B_{i,c}: M_i = sum B_{i,c}
%           1 (default) to return B_{i,c}: M_i = sum B_{i,c} * p_c
%           where
%             M_i is the class independent moment of order i
%             B_{i,c} is the class-c backward moment of order i
%             p_c is the probability of arrivals of class c
% Output:
% - MOMENTS: the backward moments as a matrix B(c,i) = B_{i,c}

% by default, moments are normalized
if nargin == 3
    NORM = 1;
end

MARKS = unique(A);

C = length(MARKS);

M = zeros(C,length(ORDERS));

for j = 1:length(ORDERS)
    k = ORDERS(j);
    for c = 1:C
        M(c,j) = mean(T.^k .* (A==MARKS(c)));
        if NORM
            M(c,j) = M(c,j) * length(T)/sum(A==MARKS(c));
        end
    end
end

end