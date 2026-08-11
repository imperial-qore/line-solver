%{
%{
 % @file pfqn_procomom.m
 % @brief ProCoMoM algorithm for computing marginal queue-length probabilities.
%}
%}

function [Pr,Q]=pfqn_procomom(L,N,Z,atol)
%{
%{
 % @brief ProCoMoM algorithm for computing marginal queue-length probabilities.
 % @fn pfqn_procomom(L, N, Z, atol)
 % @param L Service demand matrix (M x R).
 % @param N Population vector (1 x R).
 % @param Z Think time vector (1 x R).
 % @param atol Tolerance.
 % @return Pr Marginal probability matrix (M x sumN+1), Pr(k,j+1) = P(n_k = j).
 % @return Q Mean queue length vector (M x 1).
%}
%}
if nargin<3 || isempty(Z)
    Z=zeros(size(N));
end
if nargin<4
    atol=1e-14;
end
[M,R]=size(L);
sumN=sum(N);

% rescale demands per class for numerical stability
% normalized probabilities are invariant to per-class demand scaling
Lmax=max(L,[],1);
Lmax(Lmax<atol)=1;
L=L./repmat(Lmax,M,1);
Z=Z./Lmax;

% build Dn basis: multichoose(R,M) with column R zeroed, sorted by nnzpos
Dn=multichoose(R,M);
Dn(:,R)=0;
Dn=sortbynnzpos(Dn);
numDn=size(Dn,1);
basisSize=numDn*M;

% solve with auto-perturbation on rank deficiency
[Pr,rankdef]=solve_all(L,Z);
if rankdef
    % try progressively larger perturbations
    rng(23000,'twister');
    Lscale=max(abs(L(:)));
    if Lscale<atol; Lscale=1; end
    for delta_exp = [-10 -8 -6 -4]
        delta=Lscale*10^delta_exp;
        rng(23000,'twister');
        Lp=L+delta*(1+rand(M,R));
        Zp=Z+delta*(1+rand(1,R));
        [Pr_try,rd]=solve_all(Lp,Zp);
        if ~rd && all(Pr_try(:)>=-1e-6) && all(abs(sum(Pr_try,2)-1)<0.01)
            Pr=Pr_try;
            break;
        end
        Pr=Pr_try;
    end
end
Q=Pr*(0:sumN)';

    function [Pr,rankdef]=solve_all(Ls_in,Zs_in)
        L_save=L; Z_save=Z;
        L=Ls_in; Z=Zs_in;
        Pr=zeros(M,sumN+1);
        rankdef=false;
        for station=1:M
            Lrot=L;
            Lrot([station M],:)=Lrot([M station],:);
            [dist,rd]=solve_station(Lrot);
            if rd
                rankdef=true;
            end
            total=sum(dist);
            if abs(total)>0
                Pr(station,:)=dist/total;
            end
        end
        L=L_save; Z=Z_save;
    end

    function [dist,rankdef]=solve_station(Ls)
        rankdef=false;
        % pk(:,j+1) holds basis coefficients for queue length n=j
        pk=zeros(basisSize,sumN+1);
        % pexact: for empty network, all stations have probability 1 at n=0
        zero_dn=zeros(1,R);
        for kk=1:M
            pk(phash(zero_dn,kk),1)=1;
        end
        Ncur=zeros(1,R);
        for r=1:R
            for Nr=1:N(r)
                Ncur(r)=Nr;
                pklast=pk;
                pk=zeros(basisSize,sumN+1);
                [Ag,Bg,DCg,DDg]=genpmatrix(Ls,Ncur,r);
                % SVD for rank-revealing decomposition
                [U,S,V]=svd(Ag,0);
                sv=diag(S);
                tol=max(size(Ag))*eps(max(sv));
                rk=sum(sv>tol);
                if rk<basisSize
                    rankdef=true;
                    % use minimum-norm solution via truncated SVD
                    Ur=U(:,1:rk);
                    Vr=V(:,1:rk);
                    Si=diag(1./sv(1:rk));
                    pB=Vr*(Si*(Ur'*Bg));
                    pDC=Vr*(Si*(Ur'*DCg));
                    pDD=Vr*(Si*(Ur'*DDg));
                    sumNcur=sum(Ncur);
                    pk(:,1)=pB*pklast(:,1);
                    for n=1:sumNcur
                        pk(:,n+1)=pB*pklast(:,n+1)+n*pDC*pk(:,n)+n*pDD*pklast(:,n);
                    end
                else
                    % full rank: use QR for speed and accuracy
                    [Qg,Rg]=qr(Ag,0);
                    QtB=Qg'*Bg;
                    QtDC=Qg'*DCg;
                    QtDD=Qg'*DDg;
                    sumNcur=sum(Ncur);
                    pk(:,1)=Rg\(QtB*pklast(:,1));
                    for n=1:sumNcur
                        rhs=QtB*pklast(:,n+1)+n*QtDC*pk(:,n)+n*QtDD*pklast(:,n);
                        pk(:,n+1)=Rg\rhs;
                    end
                end
                % rescale to prevent overflow
                smax=max(abs(pk(:)));
                if smax>0 && isfinite(smax)
                    pk=pk/smax;
                end
            end
        end
        % extract result: first basis element at zero Dn entry
        dist=pk(phash(zero_dn,1),:);
    end

    function [A,B,DC,DD]=genpmatrix(Ls,Ncur,r)
        numRows=countrows(r);
        A=zeros(numRows,basisSize);
        B=zeros(numRows,basisSize);
        DC=zeros(numRows,basisSize);
        DD=zeros(numRows,basisSize);
        row=0;
        for d=1:numDn
            if r<=R-1 && sum(Dn(d,r:R-1))>0
                % Branch A: propagation through class boundaries
                for k=1:M
                    row=row+1;
                    colA=phash(Dn(d,:),k);
                    A(row,colA)=1;
                    if r+1<=R-1 && sum(Dn(d,r+1:R-1))>0
                        B(row,colA)=1;
                    else
                        shifted=Dn(d,:);
                        shifted(r)=shifted(r)-1;
                        colB=phash(shifted,k);
                        if colB>0 && colB<=basisSize
                            B(row,colB)=1;
                        end
                    end
                end
            else
                % Branch B: CE, PC, and extra PC equations
                if sum(Dn(d,1:r))<M
                    % CE equations: k=1..M-1
                    for k=1:M-1
                        row=row+1;
                        A(row,phash(Dn(d,:),k+1))=1;
                        A(row,phash(Dn(d,:),1))=-1;
                        for s=1:r-1
                            shifted=Dn(d,:);
                            shifted(s)=shifted(s)+1; % UP shift
                            col=phash(shifted,k+1);
                            if col>0 && col<=basisSize
                                A(row,col)=-Ls(k,s);
                            end
                        end
                        B(row,phash(Dn(d,:),k+1))=Ls(k,r);
                    end
                    % PC equations: s=1..r-1
                    for s=1:r-1
                        row=row+1;
                        nd_s=Ncur(s)-Dn(d,s);
                        A(row,phash(Dn(d,:),1))=nd_s;
                        shifted=Dn(d,:);
                        shifted(s)=shifted(s)+1; % UP shift
                        colBase=phash(shifted,1);
                        if colBase>0 && colBase<=basisSize
                            A(row,colBase)=-Z(s);
                            DC(row,colBase)=Ls(M,s);
                        end
                        for k=1:M-1
                            col=phash(shifted,k+1);
                            if col>0 && col<=basisSize
                                A(row,col)=-Ls(k,s);
                            end
                        end
                    end
                end
                % Extra PC for class r (always in Branch B)
                row=row+1;
                nd_r=Ncur(r)-Dn(d,r);
                A(row,phash(Dn(d,:),1))=nd_r;
                B(row,phash(Dn(d,:),1))=Z(r);
                for k=1:M-1
                    B(row,phash(Dn(d,:),k+1))=Ls(k,r);
                end
                DD(row,phash(Dn(d,:),1))=Ls(M,r);
            end
        end
    end

    function numRows=countrows(r)
        numRows=0;
        for d=1:numDn
            if r<=R-1 && sum(Dn(d,r:R-1))>0
                numRows=numRows+M;
            else
                if sum(Dn(d,1:r))<M
                    numRows=numRows+M+r-1;
                else
                    numRows=numRows+1;
                end
            end
        end
    end

    function col=phash(dn,i)
        pos=matchrow(Dn,dn);
        if pos<0
            col=-1;
            return;
        end
        col=(pos-1)*M+i;
    end

    function I=sortbynnzpos(I)
        for ii=1:size(I,1)-1
            for jj=ii+1:size(I,1)
                if nnzcmp(I(ii,:),I(jj,:))==1
                    v=I(ii,:);
                    I(ii,:)=I(jj,:);
                    I(jj,:)=v;
                end
            end
        end
    end

    function r=nnzcmp(i1,i2)
        nnz1=nnz(i1);
        nnz2=nnz(i2);
        if nnz1>nnz2
            r=1;
        elseif nnz1<nnz2
            r=0;
        else
            for jj=1:length(i1)
                if i1(jj)==0 && i2(jj)>0
                    r=1;
                    return
                elseif i1(jj)>0 && i2(jj)==0
                    r=0;
                    return
                end
            end
            r=0;
        end
    end

end
