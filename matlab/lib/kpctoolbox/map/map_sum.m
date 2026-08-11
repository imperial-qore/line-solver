function MAP=map_sum(MAP,n)
order = length(MAP{1});
D0=zeros(n*order);
D1=D0;

curpos=0;
for i=1:n
    D0((curpos+1):(curpos+order),(curpos+1):(curpos+order)) = MAP{1};
    if i<n
        D0((curpos+1):(curpos+order),(curpos+order+1):(curpos+2*order)) = MAP{2};
    else
        D1((curpos+1):(curpos+order),1:order) = MAP{2};
    end
    curpos = curpos + order;
end
MAP={D0,D1};
end