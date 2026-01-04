function MMAPr=mmap_timereverse(MMAP)
piq=map_pi(MMAP);
D=diag(piq);
for k=1:length(MMAP)
    MMAPr{k}=inv(D)*MMAP{k}'*D;
end
end