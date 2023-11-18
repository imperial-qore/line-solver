function MMAP = mmap_max(MMAPa, MMAPb, k)
% k is the length of the synchronisation queue

D0a = MMAPa{1}; D0b = MMAPb{1};
D1a = MMAPa{2}; D1b = MMAPb{2};

na = length(D0a); nb = length(D0b);

Ia = eye(na); Ib = eye(nb);

Za = zeros(na); Zb = zeros(nb); Z = kron(Za,Zb);

A0B0 = krons(D0a,D0b); %
A1IB = kron(D1a,Ib); %
IAB1 = kron(Ia,D1b); %
IAB0 = kron(Ia,D0b); %
A0IB = kron(D0a,Ib); %

IAB1 = kron(Ia,D1b);
A1IB = kron(D1a,Ib);

M0 = zeros(size(A0B0)*(1+k*2));
M1 = zeros(size(M0));

[iRows,iCols] = size(A0B0);

M0(1:iRows, 1:iCols*3) = [A0B0 A1IB IAB1];

for i=2:1+(k-1)*2
    r = iRows*(i-1);
    c = iCols*(i-1);
    M0(r+1:r+iRows, c+1:c+iCols) = A0B0;
end

r = iRows*(2*k-1);
c = iCols*(2*k-1);
M0(r+1:r+iRows, c+1:c+iCols) = IAB0;
r = r+iRows;
c = c+iCols;
M0(r+1:r+iRows, c+1:c+iCols) = A0IB;


for i=2:k
    r = iRows*(1+2*(i-2));
    c = iCols*(3+2*(i-2));
    M0(r+1:r+iRows, c+1:c+iCols) = A1IB;
    M0(r+1+iRows:r+iRows*2, c+1+iCols:c+iCols*2) = IAB1;
end

M1(iRows+1:iRows*3, 1:iCols) = [IAB1; A1IB];

for i=2:k
    r = iRows*(1+2*(i-1));
    c = iCols*(1+2*(i-2));
    M1(r+1:r+iRows, c+1:c+iCols) = IAB1;
    M1(r+1+iRows:r+iRows*2, c+1+iCols:c+iCols*2) = A1IB;
end

MMAP=cell(1,size(MMAPa,2));

MMAP{1} = M0;
MMAP{2} = M1;

for cls = 3:size(MMAPa,2)
    
    Mc = zeros(size(M0));
    
    D1cb = MMAPb{cls};
    D1ca = MMAPa{cls};
    
    IAB1 = kron(Ia,D1cb);
    A1IB = kron(D1ca,Ib);
    
    Mc(iRows+1:iRows*3, 1:iCols) = [IAB1; A1IB];
    
    for i=2:k
        r = iRows*(1+2*(i-1));
        c = iCols*(1+2*(i-2));
        Mc(r+1:r+iRows, c+1:c+iCols) = IAB1;
        Mc(r+1+iRows:r+iRows*2, c+1+iCols:c+iCols*2) = A1IB;
    end    
    MMAP{cls} = Mc;
end

end