function pc = mmap_pc(MMAP)
% Computes the arrival probabilities of each class for the given MMAP.
% Supports both numeric and symbolic MMAPs.

m = length(MMAP)-2;

if map_issym(MMAP)
    pc = sym(zeros(m,1));
else
    pc = zeros(m,1);
end

for i = 1:m
    pc(i) = sum(map_pie(MMAP) * (-MMAP{1}\MMAP{2+i}));
end

end