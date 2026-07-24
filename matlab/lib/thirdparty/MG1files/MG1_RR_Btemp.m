% This script multiplies B with temp, i.e., B*temp without
% storing B (only b, r1, r2, c1 and c2).
% recall B = L(b) + L(c1)L(Zr1)' + L(c2)L(Zr2)'

hd=size(temp,2);

% we start with L(Zr1)'*temp
temp1b=r1(1:(N-1)*m,:)';
temp1((N-1)*m+1:N*m,1:hd)=zeros(m,hd);
for lus1=1:N-1
  temp1((lus1-1)*m+1:lus1*m,1:hd)=temp1b(:,1:(N-lus1)*m)*temp(lus1*m+1:end,1:hd);
end
% we now obtain L(c1)*(L(Zr1)'*temp) 
for lus1=1:N
  temp1b(1:m,(lus1-1)*m+1:lus1*m)=c1((N-lus1)*m+1:(N-lus1+1)*m,:);
end    
for lus1=1:N
  temp2((lus1-1)*m+1:lus1*m,1:hd)=temp1b(:,(N-lus1)*m+1:end)*...
       temp1(1:lus1*m,:);  
end
clear temp1;
clear temp1b;
% we continue with L(Zr2)'*temp
temp1b=r2(1:(N-1)*m,:)';
temp1((N-1)*m+1:N*m,1:hd)=zeros(m,hd);
for lus1=1:N-1
  temp1((lus1-1)*m+1:lus1*m,1:hd)=temp1b(:,1:(N-lus1)*m)*temp(lus1*m+1:end,1:hd);
end
% next we store (L(c1)*L(Zr1)'+L(c2)*L(Zr2)')*temp in temp2.
for lus1=1:N
  temp1b(1:m,(lus1-1)*m+1:lus1*m)=c2((N-lus1)*m+1:(N-lus1+1)*m,:);
end    
for lus1=1:N
  temp2((lus1-1)*m+1:lus1*m,1:hd)=temp2((lus1-1)*m+1:lus1*m,1:hd)+...
      temp1b(:,(N-lus1)*m+1:end)*temp1(1:lus1*m,:);  
end
% finally we add L(b)*temp to temp2
for lus1=1:N
  temp1b(1:m,(lus1-1)*m+1:lus1*m)=b((N-lus1)*m+1:(N-lus1+1)*m,:);
end    
for lus1=1:N
  temp2((lus1-1)*m+1:lus1*m,1:hd)=temp2((lus1-1)*m+1:lus1*m,1:hd)+...
     temp1b(:,(N-lus1)*m+1:end)*temp(1:lus1*m,:); 
end
temp=temp2;
clear temp1;
clear temp2;
clear temp1b;
