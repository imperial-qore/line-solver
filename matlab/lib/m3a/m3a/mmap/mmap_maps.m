function MAPs = mmap_maps(MMAP)
% Returns K MAPs, one for each class of the MMAP[K] process.
K = length(MMAP)-2;
MAPs = cell(1,K);
for k=1:K
    MAPs{k} = {MMAP{1}+MMAP{2}-MMAP{2+k},MMAP{2+k}};
end
end