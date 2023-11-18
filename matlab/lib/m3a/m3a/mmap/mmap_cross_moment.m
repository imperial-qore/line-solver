function MC = mmap_cross_moment(mmap, k)
% Computes the k-th order moment of the inter-arrival time between an event
% of class i and an event of class j, for all possible pairs of classes.
% Input
%   mmap: the MMAP
%   k: order of the moment
% Output
%   MC: the element in (i,j) is the k-th order moment of the inter-arrival
%       time between an event of class i and an event of class j


C = length(mmap)-2;

if map_issym(mmap)
    TG = sym(zeros(C,1));
    MC = sym(zeros(C,C));
else
    TG = zeros(C,1);
    MC = zeros(C,C);
end

for i = 1:C
    if map_issym(mmap)
        TG(i) = simplify(sum(map_pie(mmap) * (-mmap{1}\mmap{2+i})));
    else
        TG(i) = sum(map_pie(mmap) * (-mmap{1}\mmap{2+i}));
    end
end

for i = 1:C
    start = map_pie(mmap) * (-mmap{1}\mmap{2+i}) / TG(i);
    for j = 1:C
        MC(i,j) = factorial(k) * sum(start * (inv(-mmap{1}))^(k+1) * mmap{2+j});
        MC(i,j) = MC(i,j) / sum(start * (-mmap{1}\mmap{2+j}));
    end
end