function E=cache_erec(gamma,m)
if nargin<3
    R=ones(length(m),length(gamma));
end
E = sub_cache_erec(gamma,m,length(gamma));
end

function E=sub_cache_erec(gamma,m,k)
h=length(m);
if sum(m)==0
    E=1;
    return
end
if sum(m)>k || min(m)<0
    E=0;
    return
end

if k==1 && sum(m)==1
    j = find(m);
    E=gamma(1,j);
    return
end
E = sub_cache_erec(gamma,m,k-1);
for j=1:h
    if m(j)>0
        E = E + gamma(k,j)*m(j)*sub_cache_erec(gamma,oner(m,j),k-1);
    end
end
end