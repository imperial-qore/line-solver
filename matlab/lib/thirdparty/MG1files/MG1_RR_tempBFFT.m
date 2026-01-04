% computes the product temp*B with B = L(b)+L(c1)L(Zr1)'+L(c2)L(Zr2)'
% temp has m rows, the result is stored in temp.
% STEP 1: temp*L(c1) in temp1a

% preparing for FFT
temp1t=reshape(temp,m^2,N)';
temp1b=c1';
temp1b=reshape(temp1b,m^2,N)';
temp1b=flipud(temp1b); % termen in omgekeerde volgorde

dim=2*2^ceil(log2(size(temp1t,1)));
temp1t=fft(temp1t,dim).';
temp1b=fft(temp1b,dim).';

% now reblock for products
temp1t=reshape(temp1t,m,m*dim);
temp1b=reshape(temp1b,m,m*dim).';

for i=1:dim
  temp1b((i-1)*m+1:i*m,1:m)=temp1t(1:m,(i-1)*m+1:i*m)*...
      temp1b((i-1)*m+1:i*m,1:m);
end
clear temp1t;
% prepare for IFFT
temp1b=temp1b.';
temp1b=reshape(temp1b,m^2,dim).';

temp1b=real(ifft(temp1b))';
temp1b=reshape(temp1b,m,m*dim)';
lushelp=N;
for lus1=1:lushelp
  temp1a(1:m,(lus1-1)*m+1:lus1*m)=temp1b((lushelp-2+lus1)*m+1:(lushelp-1+lus1)*m,:);
end  
clear temp1b;

% STEP 2: temp*L(c1)*L(Zr1)' in temp1a, thus (L(Zr1)*temp1a')'

temp1a=temp1a';

% preparing for FFT
temp1t=temp1a';
temp1t=reshape(temp1t,m^2,N)';
temp1b=[zeros(m,m); r1(1:(N-1)*m,:)]';
temp1b=reshape(temp1b,m^2,N)';

dim=2*2^ceil(log2(size(temp1t,1)));
temp1t=fft(temp1t,dim).';
temp1b=fft(temp1b,dim).';

% now reblock for products
temp1t=reshape(temp1t,m,m*dim).';
temp1b=reshape(temp1b,m,m*dim).';

for i=1:dim
  temp1b((i-1)*m+1:i*m,1:m)=temp1b((i-1)*m+1:i*m,1:m)*...
      temp1t((i-1)*m+1:i*m,1:m);
end
clear temp1t;
% prepare for IFFT
temp1b=temp1b.';
temp1b=reshape(temp1b,m^2,dim).';

temp1b=real(ifft(temp1b))';
temp1b=reshape(temp1b,m,m*dim)';

temp1a=temp1b(1:N*m,:)';
clear temp1b;

% STEP 3: temp*L(c2) in step2a

% preparing for FFT
temp1t=reshape(temp,m^2,N)';
temp1b=c2';
temp1b=reshape(temp1b,m^2,N)';
temp1b=flipud(temp1b); % termen in omgekeerde volgorde

dim=2*2^ceil(log2(size(temp1t,1)));
temp1t=fft(temp1t,dim).';
temp1b=fft(temp1b,dim).';

% now reblock for products
temp1t=reshape(temp1t,m,m*dim);
temp1b=reshape(temp1b,m,m*dim).';

for i=1:dim
  temp1b((i-1)*m+1:i*m,1:m)=temp1t(1:m,(i-1)*m+1:i*m)*...
      temp1b((i-1)*m+1:i*m,1:m);
end
clear temp1t;
% prepare for IFFT
temp1b=temp1b.';
temp1b=reshape(temp1b,m^2,dim).';

temp1b=real(ifft(temp1b))';
temp1b=reshape(temp1b,m,m*dim)';
lushelp=N;
for lus1=1:lushelp
  temp2a(1:m,(lus1-1)*m+1:lus1*m)=temp1b((lushelp-2+lus1)*m+1:(lushelp-1+lus1)*m,:);
end  
clear temp1b;
% STEP4: temp*(L(c2)*L(Zr2)'+L(c1)*L(Zr1)') in temp1a

temp2a=temp2a';

% preparing for FFT
temp1t=temp2a';
clear temp2a;
temp1t=reshape(temp1t,m^2,N)';
temp1b=[zeros(m,m); r2(1:(N-1)*m,:)]';
temp1b=reshape(temp1b,m^2,N)';

dim=2*2^ceil(log2(size(temp1t,1)));
temp1t=fft(temp1t,dim).';
temp1b=fft(temp1b,dim).';

% now reblock for products
temp1t=reshape(temp1t,m,m*dim).';
temp1b=reshape(temp1b,m,m*dim).';

for i=1:dim
  temp1b((i-1)*m+1:i*m,1:m)=temp1b((i-1)*m+1:i*m,1:m)*...
      temp1t((i-1)*m+1:i*m,1:m);
end
clear temp1t;
% prepare for IFFT
temp1b=temp1b.';
temp1b=reshape(temp1b,m^2,dim).';

temp1b=real(ifft(temp1b))';
temp1b=reshape(temp1b,m,m*dim)';

temp1a=temp1b(1:N*m,:)'+temp1a; % (temp*L(c2)*L(Zr2))' in temp1b
clear temp1b;

% STEP 5: temp*B in temp

% preparing for FFT
temp1t=reshape(temp,m^2,N)';
temp1b=b';
temp1b=reshape(temp1b,m^2,N)';
temp1b=flipud(temp1b); % termen in omgekeerde volgorde

dim=2*2^ceil(log2(size(temp1t,1)));
temp1t=fft(temp1t,dim).';
temp1b=fft(temp1b,dim).';

% now reblock for products
temp1t=reshape(temp1t,m,m*dim);
temp1b=reshape(temp1b,m,m*dim).';

for i=1:dim
  temp1b((i-1)*m+1:i*m,1:m)=temp1t(1:m,(i-1)*m+1:i*m)*...
      temp1b((i-1)*m+1:i*m,1:m);
end
clear temp1t;
% prepare for IFFT
temp1b=temp1b.';
temp1b=reshape(temp1b,m^2,dim).';

temp1b=real(ifft(temp1b))';
temp1b=reshape(temp1b,m,m*dim)';
lushelp=N;
for lus1=1:lushelp
  temp(1:m,(lus1-1)*m+1:lus1*m)=temp1b((lushelp-2+lus1)*m+1:(lushelp-1+lus1)*m,:);
end  
clear temp1b;
temp=temp+temp1a;
clear temp1a
