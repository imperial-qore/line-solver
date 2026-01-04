function [M] = mtrace_forward_moment(T,A,ORDERS,NORM)
% Computes the forward moments of a marked trace.
% Input:
% - T:      the inter-arrival times
% - C:      the class labels
% - ORDERS: vector with the orders of the moments to compute
% - NORM:   0 to return F_{i,c}: M_i = sum F_{i,c}
%           1 (default) to return F_{i,c}: M_i = sum F_{i,c} * p_c
%           where
%             M_i is the class independent moment of order i
%             F_{i,c} is the class-c forward moment of order i
%             p_c is the probability of arrivals of class c
% Output:
% - MOMENTS: the forward moments as a matrix F(c,i) = F_{i,c}

% by default, moments are normalized
if nargin == 3
    NORM = 3;
end

MARKS = unique(A);

C = length(MARKS);

M = zeros(C,length(ORDERS));

for j = 1:length(ORDERS)
    k = ORDERS(j);
    for c = 1:C
        M(c,j) = mean(T(2:end).^k .* (A(1:(end-1)) == MARKS(c)));
        if NORM
            M(c,j) = M(c,j) * length(T-1)/sum(A(1:(end-1))==MARKS(c));
        end
    end
end

end