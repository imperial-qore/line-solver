function SMMAP=mmap_sum(MMAP,n)
K=length(MMAP)-2;
for i=1:n
    MMAPs{i}=MMAP;
end
ns=[];
for i=1:n
    ns(i)=length(MMAPs{i}{1});
end
D0=zeros(sum(ns));
D1 = cell(1,K-2);
for k=2:(K+2)
    D1{k-1} = 0*D0;
end

curpos=0;
for i=1:n
    D0((curpos+1):(curpos+ns(i)),(curpos+1):(curpos+ns(i))) = MMAPs{i}{1};
    if i<n
        D0((curpos+1):(curpos+ns(i)),(curpos+ns(i)+1):(curpos+ns(i)+ns(i+1))) = MMAPs{i}{2};
    else
        for k=2:(K+2)
            D1{k-1}((curpos+1):(curpos+ns(i)),1:ns(1)) = MMAPs{i}{k};
        end
    end
    curpos = curpos + ns(i);
end
SMMAP={D0,D1{:}};
end