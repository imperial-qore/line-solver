% This script multiplies temp by B (i.e., temp*B) without
% storing B (only b, r1, r2, c1 and c2).
% recall B = L(b) + L(c1)L(Zr1)' + L(c2)L(Zr2)'

hd=size(temp,1);

% we start with temp*L(c1)
for lus1=1:N
  temp1(1:hd,(lus1-1)*m+1:lus1*m)=temp(1:hd,(lus1-1)*m+1:end)*c1(1:(N-lus1+1)*m,:);
end
% we now obtain (temp*L(c1))*L(Zr1)' 
for lus1=1:N-1
  temp1b((lus1-1)*m+1:lus1*m,1:m)=r1((N-lus1-1)*m+1:(N-lus1)*m,:)';
end    
temp2(1:hd,1:m)=zeros(hd,m);
for lus1=2:N
  temp2(1:hd,(lus1-1)*m+1:lus1*m)=temp1(1:hd,1:(lus1-1)*m)*...
       temp1b((N-lus1)*m+1:end,:);  
end
clear temp1;
clear temp1b;
% we continue with temp*Lc2
for lus1=1:N
  temp1(1:hd,(lus1-1)*m+1:lus1*m)=temp(1:hd,(lus1-1)*m+1:end)*c2(1:(N-lus1+1)*m,:);
end
% next we store temp*(L(c1)*L(Zr1)'+L(c2)*L(Zr2)') in temp2.
for lus1=1:N-1
  temp1b((lus1-1)*m+1:lus1*m,1:m)=r2((N-lus1-1)*m+1:(N-lus1)*m,:)';
end
for lus1=2:N
  temp2(1:hd,(lus1-1)*m+1:lus1*m)=temp2(1:hd,(lus1-1)*m+1:lus1*m)+...
      temp1(1:hd,1:(lus1-1)*m)*temp1b((N-lus1)*m+1:end,:);  
end
% finally we add temp*L(b) to temp2
for lus1=1:N
  temp2(1:hd,(lus1-1)*m+1:lus1*m)=temp2(1:hd,(lus1-1)*m+1:lus1*m)+...
      temp(1:hd,(lus1-1)*m+1:end)*b(1:(N-lus1+1)*m,:);
end
temp=temp2;
clear temp1;
clear temp2;
clear temp1b;