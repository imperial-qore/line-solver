% H. Emre Kankaya and Nail Akar:
% Solving Multi-Regime Feedback Fluid Queues
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% R and Q are cells with kth element (k=1...K) Q,R parameters describing the
% behaviour in regime k, Qt, Rt are cells with kth elements (k=0...K) 
% describing the behaviour at the kth boundary, T is the vector of
% thresholds of regimes (k=1...K).
%
% IMPORTANT:
% ----------
% if R{k}(m)==0 then Rt{k}(m)=0 and Rt{k+1}(m)=0 must hold
% repulsive states: boundary rate must be 0 or right or left continuous
% absorbing states: boundary rate must be 0
% all other states: boundary rate must be right or left continuous
% in all regimes the mean drift must be non-zero
%
% momnum: number of buffer length moments to compute
% arrivals: cell of input rates in all regimes. If given, buffer length
% moments will be embedded to arrival instants
%
function [pdf, pdfd, cdf, cdfm] = multiregime (Q, R, Qt, Rt, T, pdfpoints, cdfpoints)

K = length(R);
N = length(R{1});
T = [0 T];

% Convenience operations
% ~~~~~~~~~~~~~~~~~~~~~~
if length(Q)==1
    for k=2:K
        Q{k} = Q{1};
    end
end

if isempty(Qt)
    Qt{1} = Q{1};
    for k=1:K
        Qt{k+1} = Q{k};
    end
elseif length(Qt)==1
    for k=2:K+1
        Qt{k} = Qt{1};
    end
end

if isempty(Rt)
    Rt{1} = R{1};
    for k=1:K
        Rt{k+1} = R{k};
    end
end

% Obtain transfer matrices for the regimes
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
ix = 1:N;
Nnz = [];
Npos = 1;
for k=1:K
    % find negative, positive and zero states in each regime
    zix = ix(R{k}==0);
    nzix = ix(R{k}~=0);
    Nn = length(nzix);
    Qnk = Q{k}(nzix,nzix) + Q{k}(nzix,zix)*inv(-Q{k}(zix,zix))*Q{k}(zix,nzix);
    A = Qnk*diag(1./R{k}(nzix));
    [Z1,D1] = schur (A);
    [Z,D] = ordschur (Z1, D1, 5*(abs(diag(D1))<1e-10) + 2*(diag(D1)<0) + 1*(diag(D1)>0));
    zeroeig = sum(abs(diag(D1))<1e-10);
    poseig = sum(diag(D(zeroeig+1:end,zeroeig+1:end))>0);
    negeig = sum(diag(D(zeroeig+1:end,zeroeig+1:end))<0);
    
%    X1 = -D(1,2:end) / D(2:end,2:end);
%    X1b = quasitriangular (D(2:end,2:end), D(1,2:end));
    X1 = sylvester (zeros(zeroeig), -D(zeroeig+1:end, zeroeig+1:end), D(1:zeroeig, zeroeig+1:end));

    
    % the two lines below should give the same result but matlab's version
    % is much more unstable numerically 
%    X2 = lyap (-D(2:1+negeig,2:1+negeig), D(1+negeig+1:end, 1+negeig+1:end), D(2:1+negeig, 1+negeig+1:end));
    X2 = sylvester (D(zeroeig+1:zeroeig+negeig,zeroeig+1:zeroeig+negeig), -D(zeroeig+negeig+1:end, zeroeig+negeig+1:end), D(zeroeig+1:zeroeig+negeig, zeroeig+negeig+1:end));
    Y = Z * [eye(zeroeig), -X1; zeros(Nn-zeroeig,zeroeig), eye(Nn-zeroeig)] * [eye(zeroeig) zeros(zeroeig,Nn-zeroeig); zeros(negeig,zeroeig) eye(negeig) -X2; zeros(poseig,zeroeig+negeig), eye(poseig)];

    iY = inv(Y);
    iY0 = iY(1:zeroeig,:);
    iYn = iY(zeroeig+1:zeroeig+negeig,:);
    iYp = iY(zeroeig+negeig+1:end,:);
    
    At = iY*A*Y;
    An{k} = At(zeroeig+1:zeroeig+negeig, zeroeig+1:zeroeig+negeig);
    Ap{k} = At(zeroeig+negeig+1:end, zeroeig+negeig+1:end);
    
    L0{k}(:,nzix) = iY0;
    L0{k}(:,zix) = iY0*Q{k}(nzix,zix)*inv(-Q{k}(zix,zix));
    Ln{k}(:,nzix) = iYn;
    Ln{k}(:,zix) = iYn*Q{k}(nzix,zix)*inv(-Q{k}(zix,zix));
    Lp{k}(:,nzix) = iYp;
    Lp{k}(:,zix) = iYp*Q{k}(nzix,zix)*inv(-Q{k}(zix,zix));
    
    Tk = T(k+1)-T(k);
    M0{k} = [L0{k}; Ln{k}; expm(-Ap{k}*Tk)*Lp{k}];
    MT{k} = [L0{k}; expm(An{k}*Tk)*Ln{k}; Lp{k}];
    if rcond(-An{k})<1e-10 || rcond(-Ap{k})<1e-10
        An{k};
    end
    Mi{k} = [Tk*L0{k}; inv(-An{k})*(eye(negeig)-expm(An{k}*Tk))*Ln{k}; inv(Ap{k})*(eye(poseig)-expm(-Ap{k}*Tk))*Lp{k}];
    Nnz = [Nnz Nn];
    Npos = [Npos; Npos(end)+Nn];
end

% Construct and solve linear set of equations for the unknows (probability 
% masses and density parameters a0 a- a+ for all regimes)
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Neqns = (K+1)*N+sum(Nnz);

M = zeros(Neqns);
d = (K+1)*N;
p = 1;

% number of equations for the masses: (K+1)*N, corresponding variables located between indices 1...d
% number of initial densities: K*sum(Nnz), starting at d+1

% eq. (12)
M(1:N,p:p+N-1) = -Qt{1};
M(d+Npos(1):d+Npos(1)+Nnz(1)-1,p:p+N-1) = M0{1}*diag(R{1});
p = p + N;
% eq. (13)
for k=1:K-1
    M(k*N+1:(k+1)*N,p:p+N-1) = -Qt{k+1};
    M(d+Npos(k+1):d+Npos(k+1)+Nnz(k+1)-1,p:p+N-1) = M0{k+1}*diag(R{k+1});
    M(d+Npos(k):d+Npos(k)+Nnz(k)-1,p:p+N-1) = -MT{k}*diag(R{k});
    p = p + N;
end
% eq. (16)
M(K*N+1:(K+1)*N,p:p+N-1) = -Qt{K+1};
M(d+Npos(K):d+Npos(K)+Nnz(K)-1,p:p+N-1) = -MT{K}*diag(R{K});
p = p + N;
% eq. (8)
for m=1:N
    if R{1}(m)>0
        M(:,p) = zeros(Neqns,1);
        M(m,p) = 1;
        p = p + 1;
    end
end
% eq. (11)
for m=1:N
    if R{K}(m)<0
        M(:,p) = zeros(Neqns,1);
        M(K*N+m,p) = 1;
        p = p + 1;
    end
end
% eq. (9)
for k=1:K-1
    for m=1:N
        if (R{k}(m)>0 && R{k+1}(m)>0) || (R{k}(m)<0 && R{k+1}(m)<0)
            M(:,p) = zeros(Neqns,1);
            M(k*N+m,p) = 1;
            p = p + 1;
        end
    end
end
% eq. (10)
for k=1:K-1
    for m=1:N
        if (R{k}(m)<0 && R{k+1}(m)>0) && Rt{k+1}(m)~=0
            M(:,p) = zeros(Neqns,1);
            M(k*N+m,p) = 1;
            p = p + 1;
        end
    end
end
% eq. (14) 
for k=1:K-1
    for m=1:N
        if R{k}(m)<0 && Rt{k+1}(m)>=0
            M(:,p) = zeros(Neqns,1);
            M(d+Npos(k):d+Npos(k)+Nnz(k)-1,p) = MT{k}(:,m);
            p = p + 1;
        end
    end
end
% eq. (15)
for k=1:K-1
    for m=1:N
        if R{k+1}(m)>0 && Rt{k+1}(m)<=0
            M(:,p) = zeros(Neqns,1);
            M(d+Npos(k+1):d+Npos(k+1)+Nnz(k+1)-1,p) = M0{k+1}(:,m);
            p = p + 1;
        end
    end
end

% normalization
% sum up probability mass
npos = 1;
M(1:d,npos) = ones(d,1);
% sum up integrals of transfer matrices
for k=1:K
    M(d+Npos(k):d+Npos(k)+Nnz(k)-1,npos) = sum(Mi{k},2);
end

rhs = zeros(1,Neqns);
rhs(npos) = 1;

% solve linear set of equations
sol = rhs / M;

% extract results of the regimes
masses = cell(1,K+1);
a0 = cell(1,K);
an = cell(1,K);
ap = cell(1,K);
masses{1} = sol(1:N);
for k=1:K
    masses{k+1} = sol(k*N+1:(k+1)*N);
    avec = sol(d+Npos(k):d+Npos(k)+Nnz(k)-1);
    a0{k} = avec(1:Nnz(k)-size(An{k},1)-size(Ap{k},1));
    an{k} = avec(length(a0{k})+1:length(a0{k})+size(An{k},1));
    ap{k} = avec(length(a0{k})+length(an{k})+1:length(a0{k})+length(an{k})+size(Ap{k},1));
end

% calculate pdf at the requested points
pdf = [];
for p=pdfpoints
    k=0;
    while k<K && p>=T(k+1)
        k = k + 1;
    end
    presp = a0{k}*L0{k} + an{k}*expm(An{k}*(p-T(k)))*Ln{k} + ap{k}*expm(-Ap{k}*(T(k+1)-p))*Lp{k};
    pdf = [pdf; presp];
end

pdfd = [];
for p=pdfpoints
    k=0;
    while k<K && p>=T(k+1)
        k = k + 1;
    end
    presp = an{k}*An{k}*expm(An{k}*(p-T(k)))*Ln{k} + ap{k}*Ap{k}*expm(-Ap{k}*(T(k+1)-p))*Lp{k};
    pdfd = [pdfd; presp];
end

% calculate cdf at the requested points
cdf = [];
cdfm = [];
for ix=1:length(cdfpoints)
    c=cdfpoints(ix);
    cres = zeros(1,N);
    cresm = zeros(1,N);
    k=0;
    while k<K && c>=T(k+1)
        if k>0
            cres = cres + [a0{k} an{k} ap{k}] * Mi{k};
            cresm = cresm + [a0{k} an{k} ap{k}] * Mi{k};
        end
        cresm = cresm + masses{k+1};
        if c>T(k+1)
            cres = cres + masses{k+1};
        end
        k = k + 1;
    end
    if k==K && c==T(k+1)
        cresm = cresm + masses{k+1};       
    end
    crem = c - T(k);
    Tk = T(k+1)-T(k);
    val = a0{k}*L0{k}*crem + an{k}*inv(-An{k})*(eye(size(An{k}))-expm(An{k}*crem))*Ln{k} + ap{k}*inv(-Ap{k})*(expm(-Ap{k}*Tk) - expm(-Ap{k}*(Tk-crem)))*Lp{k};
    cres = cres + val;
    cresm = cresm + val;
    cdf = [cdf; cres];
    cdfm = [cdfm; cresm];
end

