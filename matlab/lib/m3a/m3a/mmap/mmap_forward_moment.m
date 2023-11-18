function MOMENTS = mmap_forward_moment(MMAP,ORDERS,NORM)
% Computes the theoretical forward moments of an MMAP.
% Input:
% - MMAP:   the MMAP
% - ORDERS: vector with the orders of the moments to compute
% - NORM:   0 to return F_{i,c}: M_i = sum F_{i,c}
%           1 (default) to return F_{i,c}: M_i = sum F_{i,c} * p_c
%           where
%             M_i is the class independent moment of order i
%             F_{i,c} is the class-c forward moment of order i
%             p_c is the probability of arrivals of class c
% Output:
% - MOMENTS: the forward moments as a matrix F(c,i) = F_{i,c}

% by default moments are normalized
if nargin == 2
    NORM = 1; 
end

C = length(MMAP)-2;
K = length(ORDERS);

if (map_issym(MMAP))
    MOMENTS = sym(zeros(C,K));
else
    MOMENTS = zeros(C,K);
end

pie = map_pie(MMAP);
M = inv(-MMAP{1});

for a = 1:C
    if NORM
        pa = sum(pie * (-MMAP{1}\MMAP{2+a}));
    else
        pa = 1;
    end
    for h = 1:length(ORDERS)
        k = ORDERS(h);
        fk = factorial(k);
		MOMENTS(a,h) = fk/pa * sum(pie * (-MMAP{1} \ MMAP{2+a}) * M^k);
    end
end

end