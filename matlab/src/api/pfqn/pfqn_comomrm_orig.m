%{
%{
 % @file pfqn_comomrm_orig.m
 % @brief Original CoMoM implementation for finite repairman model.
%}
%}

%{
%{
 % @brief Original CoMoM implementation for finite repairman model.
 % @fn pfqn_comomrm_orig(L, N, Z, atol)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param atol Absolute tolerance for numerical computations.
 % @return lG Logarithm of normalizing constant.
%}
%}
function lG=pfqn_comomrm_orig(L,N,Z,atol)
% comom for a finite repairment model
if size(L,1)~=1
    line_error(mfilename,'The solver accepts at most a single queueing station.')
end
if nargin<4
    m=1;
end
lambda = 0*N;
[~,L,N,Z,lG0] = pfqn_nc_sanitize(lambda,L,N,Z,atol);
[M,R]=size(L);
% rescale demands
Lmax = L; % use L
Lmax(Lmax<atol)=Z(Lmax<atol); % unless zero
L = L./repmat(Lmax,M,1);
Z = Z./repmat(Lmax,M,1);
% sort from smallest to largest
[~,rsort] = sort(Z,'ascend');
L=L(:,rsort);
Z=Z(:,rsort);
% initialize
nvec=zeros(1,R);
h=ones(2,1);
lh=log(h) + factln(sum(nvec)+M-1) - sum(factln(nvec));
h=exp(lh);
scale=zeros(1,sum(N));
% iterate
for r=1:R
    for Nr=1:N(r)
        nvec(r)=nvec(r)+1;
        if Nr==1
            if r>1
                P = zeros(2*(r-1),2*r);
                r1 = r-1;
                P(1:r1,1:r1) = eye(r1);
                P((r1+1):2*r1,(r+1):(2*r-1)) = eye(r1);
                h1=h'*P;
                h1=h1';
                h1(r)=h_1(1)*nvec(r1)/(sum(nvec)-1)/scale(nt);
                h1(end)=h_1(r1+1)*nvec(r1)/(sum(nvec)-1)/scale(nt);
                h=h1;
            end
            A = zeros(2*r);
            DA = zeros(2*r);
            B = zeros(2*r);
            % 1 CE for G+
            A(1,1) = 1;
            A(1,2:r) = -L(1,1:r-1);
            A(1,r+1) = -1;
            B(1,1) = L(1,r);
            % Class-1..(R-1) PCs for G
            for s=1:(r-1)
                A(1+s,r+1) = N(s);
                A(1+s,r+1+s) = -Z(s);
                A(1+s,1+s) = -m*L(1,s);
            end
            % Class-R PCs for G and Gr (r=1...R-1)
            A(r+1:2*r,r+1:2*r) = Nr*eye(r);
            DA(r+1:2*r,r+1:2*r) = eye(r);
            B(r+1:2*r,1:r) = m*L(1,r)*eye(r);
            B(r+1:2*r,r+1:2*r) = Z(r)*eye(r);
            C = A(1:r,1:r);
            A12 = A(1:r,r+1:2*r);
            B1r = B(1:r,:);
            B2r = B(r+1:2*r,:);
            F1r = [inv(C)*B1r; 0*B2r];
            F2r = [-C\A12*B2r; B2r];
        end
        h_1 = h;
        h = (nvec(r)*F1r+F2r)*h_1/(sum(nvec)+M-1);
        nt=sum(nvec);
        scale(nt)=abs(sum(sort(h)));
        h = abs(h)/scale(nt); % rescale so that |h|=1
    end
end
% unscale and return the log of the normalizing constant
lG=lG0+log(h(end-(R-1))) + factln(sum(N)+M-1) - sum(factln(N)) + N*log(Lmax)' + sum(log(scale));
end
