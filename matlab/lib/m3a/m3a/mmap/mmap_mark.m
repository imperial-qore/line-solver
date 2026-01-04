function mmap = mmap_mark(MMAP, prob)
% MMAP_MARK. Marks arrivals from a MMAP according to give probabilities.
% MMAP = MMAP_MARK(MMAP, PROB) takes a MMAP with K types and a probability
% matrix with element PROB(k,s) giving the probability that a type-k
% arrival should be marked as a class-s arrival and returns a new MMAP
% that has R classes as output types.
[K,R]  = size(prob);
mmap = cell(1,2+R);
mmap{1} = MMAP{1};
mmap{2} = MMAP{2};
for r=1:R
    mmap{2+r} = zeros(length(MMAP{1}));
    for k=1:K
        mmap{2+r} = mmap{2+r} + MMAP{2+k}*prob(k,r);
    end
end
end