function idc = mmap_count_idc(mmap, t)
% Computes the per-class Index of Dispersion of Counts for the given MMAP
% at resolution t.

m = mmap_count_mean(mmap,t);
v = mmap_count_var(mmap,t);
idc = v ./ m;

end