function IDC = mmap_idc(MMAP)
tinf = 1e6 / sum(mmap_lambda(MMAP));
IDC = mmap_count_var(MMAP,tinf)/mmap_count_mean(MMAP,tinf);
end