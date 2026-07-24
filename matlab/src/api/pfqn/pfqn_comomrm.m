%{
%{
 % @file pfqn_comomrm.m
 % @brief CoMoM (Class-Oriented Method of Moments) for finite repairman model.
%}
%}

%{
%{
 % @brief CoMoM (Class-Oriented Method of Moments) for finite repairman model.
 % @fn pfqn_comomrm(L, N, Z, m, atol)
 % @param L Service demand matrix.
 % @param N Population vector.
 % @param Z Think time vector.
 % @param m Replication factor (default: 1).
 % @param atol Absolute tolerance for numerical computations.
 % @return lG Logarithm of normalizing constant.
 % @return lGbasis Logarithm of basis functions.
%}
%}
function [lG,lGbasis]=pfqn_comomrm(L,N,Z,m,atol)
% comom for a finite repairment model
[M,R]=size(L);
if M~=1
    line_error(mfilename,'The solver accepts at most a single queueing station.')
end
if nargin<4
    m=1;
end
% Same omission as pfqn_comom: atol was read by pfqn_nc_sanitize below but
% never defaulted, so the 3- and 4-argument forms in the header were dead.
if nargin<5
    atol=1e-14;
end
lambda = 0*N;
[~,L,N,Z,lG0] = pfqn_nc_sanitize(lambda,L,N,Z,atol);
% R must be re-read: sanitize DROPS zero-population and zero-demand classes,
% so the pre-sanitize R overruns the shortened Z, L and N below.
R = size(L,2);
zerothinktimes = zeros(1,R);
numZeroThinkTimes = 0;
for r = 1:R
    if Z(r) < GlobalConstants.FineTol
        % see _kb/03-api-layer.md (pfqn/ family: scaling, log-domain switches, dispatch gates)
        numZeroThinkTimes = numZeroThinkTimes + 1;
        zerothinktimes(numZeroThinkTimes) = r;
    end
end
% initialize
nvec=zeros(1,R);
if numZeroThinkTimes > 0
    for z = 1:numZeroThinkTimes
        idx = zerothinktimes(z);
        nvec(idx) = N(idx);
    end
    lh=zeros(2+2*numZeroThinkTimes,1);
    lhIdx = 1;
    % these are trivial models with a single queueing station with demands all equal to one and think time 0
    lh(lhIdx,1) = (factln(sum(nvec)+m+1-1)-sum(factln(nvec)));
    lhIdx = lhIdx + 1;
    for z = 1:numZeroThinkTimes
        s = zerothinktimes(z);
        nvec_s = oner(nvec,s);
        lh(lhIdx,1) = (factln(sum(nvec_s)+m+1-1)-sum(factln(nvec_s)));
        lhIdx = lhIdx + 1;
    end
    lh(lhIdx,1) = (factln(sum(nvec)+m-1)-sum(factln(nvec)));
    lhIdx = lhIdx + 1;
    for z = 1:numZeroThinkTimes
        s = zerothinktimes(z);
        nvec_s = oner(nvec,s);
        lh(lhIdx,1) = (factln(sum(nvec_s)+m-1)-sum(factln(nvec_s)));
        lhIdx = lhIdx + 1;
    end
else
    lh=zeros(2,1);
end
h=exp(lh);
if numZeroThinkTimes==R
    lGbasis = log(h);
    lG = lG0 + log(h(end-R));
    return
else
    scale = ones(1,sum(N));
    nt = sum(nvec);
    h_1=h;
    %iterate
    for r=(numZeroThinkTimes+1):R
        F1r = zeros(2*r);
        F2r = zeros(2*r);
        for Nr=1:N(r)
            nvec(r)=nvec(r)+1;
            if Nr==1
                if r> numZeroThinkTimes+1
                    hr = zeros(2*r,1);
                    hr(1:(r-1)) = h(1:(r-1));
                    hr((r+1):(2*r-1)) = h(((r-1)+1):2*(r-1));
                    h=hr;
                    % update scalings
                    if nt>0
                        h(r)=h_1(1)/scale(nt);
                        h(end)=h_1((r-1)+1)/scale(nt);
                    end
                end
                % CE for G+
                A12 = zeros(r);
                A12(1,1) = -1;
                % Class-1..(R-1) PCs for G
                for s=1:(r-1)
                    A12(1+s,1) = N(s);
                    A12(1+s,1+s) = -Z(s);
                end
                % Class-R PCs
                B2r = [m*L(1,r)*eye(r), Z(r)*eye(r)];
                % explicit formula for inv(C)
                iC=-eye(r)/m;
                iC(1,:)=-1/m;
                iC(1)=1;
                % explicit formula for F1r
                F1r = zeros(2*r); F1r(1,1)=1;
                % F2r by the definition
                F2r = [-iC*A12*B2r; B2r];
            end
            h_1 = h;
            h = (F1r+F2r/nvec(r))*h_1;
            nt = sum(nvec);
            scale(nt) = abs(sum(sort(h)));
            h = abs(h)/scale(nt); % rescale so that |h|=1
        end
    end

    % unscale and return the log of the normalizing constant
    lG = lG0 + log(h(end-(R-1))) + sum(log(scale));
    lGbasis = log(h)  + sum(log(scale));
end
end
