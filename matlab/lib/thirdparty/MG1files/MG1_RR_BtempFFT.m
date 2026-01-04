% computes product B*temp with B = L(b)+L(c1)L(Zr1)'+L(c2)L(Zr2)'
% temp has m columns, the result is stored in temp.
% STEP 1: L(Zr1)'*temp in temp1a, thus (temp'*L(Zr1))' 

temp=temp';

% preparing for FFT
temp1t=reshape(temp,m^2,N)';
temp1b=[zeros(m,m); r1(1:m*(N-1),:)]'; % Zr1'
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
% prepare for IFFT
temp1b=temp1b.';
temp1b=reshape(temp1b,m^2,dim).';

temp1b=real(ifft(temp1b))';
temp1b=reshape(temp1b,m,m*dim)';
lushelp=N;
for lus1=1:lushelp
  temp1a(1:m,(lus1-1)*m+1:lus1*m)=temp1b((lushelp-2+lus1)*m+1:(lushelp-1+lus1)*m,:);
end  

temp1a=temp1a';

% STEP 2: L(c1)L(Zr1)'temp in temp1a

% preparing for FFT
temp1t=temp1a';
temp1t=reshape(temp1t,m^2,N)';
temp1b=c1';
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

temp1a=temp1b(1:m*N,:);
clear temp1b;

% STEP 3: LZr2'*temp in temp2a

% preparing for FFT
temp1t=reshape(temp,m^2,N)';
temp1b=[zeros(m,m); r2(1:m*(N-1),:)]'; % Zr1'
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
% prepare for IFFT
temp1b=temp1b.';
temp1b=reshape(temp1b,m^2,dim).';

temp1b=real(ifft(temp1b))';
temp1b=reshape(temp1b,m,m*dim)';
lushelp=N;
for lus1=1:lushelp
  temp2a(1:m,(lus1-1)*m+1:lus1*m)=temp1b((lushelp-2+lus1)*m+1:(lushelp-1+lus1)*m,:);
end  

temp2a=temp2a';

% STEP 4: (L(c1)*L(Zr1)'+L(c2)*L(Zr2)')*temp in temp1a

% preparing for FFT
temp1t=temp2a';
clear temp2a;
temp1t=reshape(temp1t,m^2,N)';
temp1b=c2';
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

temp1a=temp1a+temp1b(1:m*N,:);
clear temp1b;

% STEP 5: B*temp in temp.

temp=temp'; % transposing of step 1 is restored

% preparing for FFT
temp1t=temp';
temp1t=reshape(temp1t,m^2,size(temp,1)/m)';
temp1b=b';
temp1b=reshape(temp1b,m^2,size(b,1)/m)';

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

temp=temp1a+temp1b(1:size(temp,1),:);
clear temp1b;
clear temp1a;
