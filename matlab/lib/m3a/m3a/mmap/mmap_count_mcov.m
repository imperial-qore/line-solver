function S = mmap_count_mcov(mmap, t)
% Computes the count covariance between each pair of classes at time scale
% t.
% INPUT
% - mmap: the marked MAP with m classes
% - t: the time scale
% OUTPUT
% - S: the m x m covariance matrix

m = size(mmap,2)-2;

% on the diagonal we have the per-class variance
mV = mmap_count_var(mmap, t);
S = diag(mV);

for i = 1:m
    for j = 1:m
        if i ~= j
            % compute variance between this pair of classes
            mmap2 = cell(1,4);
            mmap2{1} = mmap{1};
            mmap2{2} = mmap{2};
            mmap2{3} = mmap{2+i} + mmap{2+j};
            mmap2{4} = mmap2{2}-mmap2{3};
            pV = mmap_count_var(mmap2, t);
            S(i,j) = 1/2*(pV(1) - mV(i) - mV(j));
        end
    end
end

end