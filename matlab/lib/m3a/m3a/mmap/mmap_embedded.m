function Pc = mmap_embedded(MMAP)
K = length(MMAP) - 2;
Pc = cell(1,K);
for k=1:K
    Pc{k} = -MMAP{1} \ MMAP{2+k};
end
end