function MAP=map_sumind(MAPs)
n=length(MAPs); % MAPs to be summed
order=[];
for i=1:n
    order(i)=length(MAPs{i}{1});
end
D0=zeros(sum(order));
D1=D0;

curpos=0;
for i=1:n
   D0((curpos+1):(curpos+order(i)),(curpos+1):(curpos+order(i))) = MAPs{i}{1};       
   if i<n
       D0((curpos+1):(curpos+order(i)),(curpos+order(i)+1):(curpos+order(i)+order(i+1))) = MAPs{i}{2}*ones(order(i),1)*map_pie(MAPs{1+i});
   else       
       D1((curpos+1):(curpos+order(i)),1:order(1)) = MAPs{i}{2}*ones(order(i),1)*map_pie(MAPs{1});
   end
   curpos = curpos + order(i);
end
MAP={D0,D1};
end