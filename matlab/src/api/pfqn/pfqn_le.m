%{
%{
 % @file pfqn_le.m
 % @brief Logistic expansion (LE) asymptotic approximation for normalizing constant.
%}
%}

%{
%{
 % @brief Logistic expansion (LE) asymptotic approximation for normalizing constant.
 % @fn pfqn_le(L, N, Z)
 % @param L Service demand matrix (MxR).
 % @param N Population vector (1xR).
 % @param Z Think time vector (1xR).
 % @return Gn Estimated normalizing constant.
 % @return lGn Logarithm of normalizing constant.
%}
%}
function [Gn,lGn]=pfqn_le(L,N,Z)
% [GN,LGN]=PFQN_LE(L,N,Z)

% PFQN_LE Asymptotic solution of closed product-form queueing networks by
% logistic expansion
%
% [Gn,lGn]=pfqn_le(L,N,Z)
% Input:
% L : MxR demand matrix. L(i,r) is the demand of class-r at queue i
% N : 1xR population vector. N(r) is the number of jobs in class r
% Z : 1xR think time vector. Z(r) is the total think time of class r
%
% Output:
% Gn : estimated normalizing constat
% lGn: logarithm of Gn. If Gn exceeds the floating-point range, only lGn
%      will be correctly estimated.
%
% Reference:
% G. Casale. Accelerating performance inference over closed systems by
% asymptotic methods. ACM SIGMETRICS 2017.
% Availble at: http://dl.acm.org/citation.cfm?id=3084445


[M,R]=size(L);

if isempty(L) || isempty(N) || sum(N)==0 || sum(L(:))<1e-4
    lGn = - sum(factln(N)) + sum(N.*log(sum(Z,1)));
    Gn=exp(lGn);
elseif nargin<3%~exist('Z','var')
    umax=pfqn_le_fpi(L,N);
    A=pfqn_le_hessian(L,N,umax'); % slightly faster than pfqn_le_hessianZ
    S=0; 
	for r=1:R 
		S=S+N(r)*log(umax'*L(:,r)); 
	end
    lGn = multinomialln([N,M-1]) + factln(M-1) + (M-1)*log(sqrt(2*pi)) - log(sqrt(det(A))) + sum(log(umax)) + S;
    Gn=exp(lGn);
else % Z>0
    [umax,vmax]=pfqn_le_fpiZ(L,N,Z);
    A=pfqn_le_hessianZ(L,N,Z,umax',vmax);
    S=0; 
	for r=1:R 
		S=S+N(r)*log(Z(r)+vmax*umax'*L(:,r)); 
	end
    lGn = -sum(factln(N)) -vmax + M*log(vmax) + M*log(sqrt(2*pi)) - log(sqrt(det(A))) + sum(log(umax)) + S;
    Gn=exp(lGn);
end
end

function [u,d]=pfqn_le_fpi(L,N)
% [U,D]=PFQN_LE_FPI(L,N)

% find location of mode of gaussian
[M,R]=size(L);
u=ones(M,1)/M;
u_1=Inf*u;
d=[];
while norm(u-u_1,1)>1e-10
    u_1=u;
    for i=1:M
        u(i)=1/(sum(N)+M);
        for r=1:R
            u(i)=u(i)+N(r)/(sum(N)+M)*L(i,r)*u_1(i)/(u_1'*L(:,r));
        end
    end
    d(end+1,:)=abs(u-u_1)';
end
end

function [u,v,d]=pfqn_le_fpiZ(L,N,Z)
% [U,V,D]=PFQN_LE_FPIZ(L,N,Z)

% find location of mode of gaussian
[M,R]=size(L);
eta = sum(N)+M;
u=ones(M,1)/M;
v=eta+1;
u_1=Inf*u;
v_1=Inf*v; %#ok<NASGU>
d=[];
while norm(u-u_1,1)>1e-10
    u_1=u;
    v_1=v;
    for ist=1:M
        u(ist)=1/eta;
        for r=1:R
            u(ist)=u(ist)+(N(r)/eta)*(Z(r)+v*L(ist,r))*u_1(ist)/(Z(r)+v*u_1'*L(:,r));
        end
    end
    for r=1:R
        xi(r)=N(r)/(Z(r)+v*u_1(:)'*L(:,r));
    end
    v=eta+1;
    for r=1:R
        v=v-xi(r)*Z(r);
    end
    d(end+1,:)=abs(u-u_1)'+abs(v-v_1);
end

end

function hu=pfqn_le_hessian(L,N,u0)
% HU=PFQN_LE_HESSIAN(L,N,U0)

% find hessian of gaussian
[M,R]=size(L);
Ntot=sum(N);
hu=zeros(M-1);
for i=1:(M-1)
    for j=1:(M-1)
        if i~=j
            hu(i,j)=-(Ntot+M)*u0(i)*u0(j);
            for r=1:R
                hu(i,j)=hu(i,j)+N(r)*L(i,r)*L(j,r)*(u0(i)*u0(j))/(u0*L(:,r))^2;
            end
        else % i=j
            hu(i,j)=(Ntot+M)*u0(i)*sum(allbut(u0,i));
            for r=1:R
                hu(i,j)=hu(i,j)-N(r)*L(i,r)*u0(i)*(allbut(u0,i)*L(allbut(1:M,i),r))/(u0*L(:,r))^2;
            end
        end
    end
end
end

function A=pfqn_le_hessianZ(L,N,Z,u,v)
% A=PFQN_LE_HESSIANZ(L,N,Z,U,V)

% find hessian of gaussian
[K,R]=size(L);
Ntot=sum(N);
A=zeros(K);
csi = zeros(1,R);
for r=1:R
    csi(r)=N(r)/(Z(r)+v*u*L(:,r));
end
Lhat = zeros(K,R);
for k=1:K
    for r=1:R
        Lhat(k,r)=Z(r)+v*L(k,r);
    end
end
eta=Ntot+K;
for i=1:K
    for j=1:K
        if i~=j
            A(i,j)=-eta*u(i)*u(j);
            for r=1:R
                A(i,j)=A(i,j)+csi(r)^2*Lhat(i,r)*Lhat(j,r)*(u(i)*u(j))/N(r);
            end
        end
    end
end
for i=1:K
    A(i,i)=-sum(allbut(A(i,:),i));
end
A=A(1:(K-1),1:(K-1));
A(K,K)=1;
for r=1:R
    A(K,K)=A(K,K)-(csi(r)^2/N(r))*Z(r)*u*L(:,r);
end
A(K,K)=v*A(K,K);
for i=1:(K-1)
    A(i,K)=0;
    for r=1:R
        A(i,K)=A(i,K)+v*u(i)*((csi(r)^2/N(r))*Lhat(i,r)*(u*L(:,r))-csi(r)*L(i,r));
    end
    A(K,i)=A(i,K);
end
end

function y=allbut(y,xset)
% Y=ALLBUT(Y,XSET)

y=y(setdiff(1:length(y),xset));
end

function mln=multinomialln(n)
% MLN=MULTINOMIALLN(N)

mln = factln(sum(n))- sum(factln(n));
end

function lf=factln(n)
% LF=FACTLN(N)

lf = gammaln(1+n);
end
