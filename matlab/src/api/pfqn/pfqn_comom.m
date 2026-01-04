%{
%{
 % @file pfqn_comom.m
 % @brief CoMoM algorithm for computing the normalizing constant.
%}
%}

function lG=pfqn_comom(L,N,Z,atol)
%{
%{
 % @brief CoMoM algorithm for computing the normalizing constant.
 % @fn pfqn_comom(L, N, Z, atol)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param atol Tolerance.
 % @return lG Logarithm of the normalizing constant.
%}
%}
if nargin<3
    Z=0*N;
end
[M,R]=size(L);
% rescale demands
Lmax = L; % use L
Lmax(Lmax<atol)=Z(Lmax<atol); % unless zero
Lmax = max(Lmax,[],1);
L = L./repmat(Lmax,M,1);
Z = Z./repmat(Lmax,M,1);
% sort from smallest to largest
%[~,rsort] = sort(Z,'ascend');
%L=L(:,rsort);
%Z=Z(:,rsort);
% prepare comom data structures
Dn=multichoose(R,M);
Dn(:,R)=0;
Dn=sortbynnzpos(Dn);
% initialize
nvec=zeros(1,R);
h=ginit(L);
lh=log(h) + factln(sum(nvec)+M-1) - sum(factln(nvec));
h=exp(lh);
scale=zeros(1,sum(N));
% iterate
for r=1:R
    for Nr=1:N(r)
        nvec(r)=nvec(r)+1;
        if Nr==1
            [A,B,DA]=genmatrix(L,nvec,Z,r);
        else
            A=A+DA;
        end
        b = B*h*nvec(r)/(sum(nvec)+M-1);
        h = A\b;
        nt=sum(nvec);
        scale(nt)=abs(sum(sort(h)));
        h = abs(h)/scale(nt); % rescale so that |h|=1
    end
end
% unscale and return the log of the normalizing constant
lG=log(h(end-(R-1))) + factln(sum(N)+M-1) - sum(factln(N)) + N*log(Lmax)' + sum(log(scale));

    function [A,B,DA,rr]=genmatrix(L,N,Z,r)
        [M,R]=size(L);
        A=zeros(nchoosek(M+R-1,M)*(M+1));
        DA=zeros(nchoosek(M+R-1,M)*(M+1));
        B=zeros(nchoosek(M+R-1,M)*(M+1));
        row=0;
        rr=[]; % artificial rows
        lastnnz=0;
        for d=1:length(Dn)
            hnnz=hashnnz(Dn(d,:),R);
            if hnnz~=lastnnz
                lastnnz=hnnz;
            end
        end
        for d=1:length(Dn)
            if sum(Dn(d,(r):R-1))>0
                % dummy rows for unused norm consts
                for k=0:M
                    row=row+1;
                    col = hash(N,N-Dn(d,:),k+1);
                    A(row,col)=1;
                    if sum(Dn(d,(r+1):R-1))>0
                        col = hash(N,N-Dn(d,:),k+1);
                        B(row,col)=1;
                    else
                        er=zeros(1,R); er(r)=1;
                        col = hash(N,N-Dn(d,:)+er,k+1);
                        B(row,col)=1;
                    end
                end
            else
                if sum(Dn(d,1:r))<M
                    for k=1:M
                        % add CE
                        row=row+1;
                        A(row,hash(N,N-Dn(d,:),k+1))=1;
                        A(row,hash(N,N-Dn(d,:),0+1))=-1;
                        for s=1:r-1
                            A(row,hash(N,oner(N-Dn(d,:),s),k+1))=-L(k,s);
                        end
                        B(row,hash(N,N-Dn(d,:),k+1))=L(k,r);
                    end
                    for s=1:(r-1)
                        % add PC to A
                        row=row+1;
                        n=N-Dn(d,:);
                        A(row,hash(N,n,0+1))=n(s);
                        A(row,hash(N,oner(n,s),0+1))=-Z(s);
                        for k=1:M
                            A(row,hash(N,oner(n,s),k+1))=-L(k,s);
                        end
                        B(row,:)=0;
                    end
                end
            end
        end
        %add PC of class R
        for d=1:length(Dn)
            if sum(Dn(d,(r):R-1))<=0
                row=row+1;
                n=N-Dn(d,:);
                A(row,hash(N,n,0+1))=n(r);
                DA(row,hash(N,n,0+1))=1;
                B(row,hash(N,n,0+1))=Z(r);
                for k=1:M
                    B(row,hash(N,n,k+1))=L(k,r);
                end
            end
        end
        rr=unique(rr);
    end

    function val=hashnnz(dn,R)
        val=0;
        for t=1:R
            if dn(t)==0
                val=val+2^(t-1);
            end
        end
    end

    function col=hash(N,n,i)
        if i==1
            col=size(Dn,1)*M+matchrow(Dn,N-n);
        else
            col=(matchrow(Dn,N-n)-1)*M+i-1;
        end
    end

    function g=ginit(L)
        e1=zeros(1,R);
        g=zeros(size(Dn,1)*(M+1),1);
        for i=0:M
            g(hash(N,N,i+1))=1;
        end
        g=g(:);
    end

    function I=sortbynnzpos(I)
        % sorts a set of combinations with repetition according to the number of
        % nonzeros
        for i=1:size(I,1)-1
            for j=i+1:size(I,1)
                if nnzcmp(I(i,:),I(j,:))==1
                    v=I(i,:);
                    I(i,:)=I(j,:);
                    I(j,:)=v;
                end
            end
        end
    end

    function r=nnzcmp(i1,i2) % return 1 if i1<i2
        nnz1=nnz(i1);
        nnz2=nnz(i2);
        if(nnz1>nnz2)
            r=1; % i2 has more zeros and is thus greater
        elseif(nnz1<nnz2)
            r=0; % i1 has more zeros and is thus greater
        else %nnz1==nnz2
            for j=1:length(i1)
                if i1(j)==0 & i2(j)>0
                    r=1; % i2 has the left-most zero
                    return
                elseif i1(j)>0 & i2(j)==0
                    r=0;
                    return
                end
            end
            r=0;
        end
    end

end