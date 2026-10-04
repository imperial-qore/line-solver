function [M] = mtrace_moment(T,A,ORDERS,AFTER,NORM)
% Computes the empirical class-dependent moments of a multi-class trace.
% T:      vector of inter-arrival times
% A:      vector of class labels
% ORDERS: vector with the orders of the moments to compute
% AFTER:  0 to compute moments of Horvath variables
%         1 to compute moments of Bucholz variables
% NORM:   0 to not normalize, M_i = sum M_{i,c}
%         1 to normalize, M_i = sum M_{i,c} * g_c
%         where M_i is the class independent moment of order i
%               M_{i,c} is the class c moment of order i
%               g_c is the fraction of arrivals of class c

if nargin == 3
    AFTER = 0;
end

if nargin <= 4
    NORM = 0;
end

MARKS = unique(A);

C = length(MARKS);

M = zeros(C,length(ORDERS));

for j = 1:length(ORDERS)
    k = ORDERS(j);
    for c = 1:C
        if AFTER
            M(c,j) = mean(T(2:end).^k .* (A(1:(end-1)) == MARKS(c)));
        else
            M(c,j) = mean(T.^k .* (A==MARKS(c)));
        end
        if NORM
            if AFTER
                M(c,j) = M(c,j) * length(T-1)/sum(A(1:(end-1))==MARKS(c));
            else
                M(c,j) = M(c,j) * length(T)/sum(A==MARKS(c));
            end
        end
    end
end

end