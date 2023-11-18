function [ql,qlT,soj_alpha,wait_alpha,Smat]=Q_CT_MMAPK_PHK_1(D0,D,alpha,S,varargin)
%[ql,qlT,soj_alpha,wait_alpha,Smat]=Q_CT_MMAPK_PHK_1(D0,D,alpha,S) 
%   computes the  Sojourn time and Waiting time distribution of 
%   a continuous-time MMAP[K]/PH[K]/1/FCFS queue (per type and overall)
%   
%   INPUT PARAMETERS:
%   * MMAP[K] arrival process (with m states)
%     D{i} holds the mxm matrix D_i, together with the mxm matrix D0 
%     these characterize the MMAP[K] arrival process
%
%   * PH[K] services time distributions (i=1:K)
%     alpha{i} holds the 1xmi alpha vector of the PH service time of
%     type i customers
%     S{i} holds the mixmi matrix S of the PH service time of a type
%     i customer
%
%   RETURN VALUES: 
%   * Per type queue length distribution, 
%     ql{k}(i) = Prob[(i-1) type k customers in the queue]
%   * Overall queue length distribution,
%     qlT(i) = Prob[(i-1) customers in the queue (any type)]
%   * Overall sojourn time distribution, 
%     is PH characterized by (soj_alpha{K+1},Smat)
%   * Type k sojourn time distribution,
%     is PH characterized by (soj_alpha{k},Smat)
%   * Overall waiting time distribution, 
%     is PH characterized by (wait_alpha{K+1},Smat)
%   * Type k waiting time distribution,
%     is PH characterized by (wait_alpha{k},Smat)
%
%   USES: Q_Sylvest
%
%   Optional Parameters:
% 
%       Mode: 'Sylves' solves a Sylvester matrix equation at each step 
%             using an Hessenberg algorithm
%             'Direct' solves the Sylvester matrix equation at each
%             step by rewriting it as a (large) system of linear equations
%             (default: 'Sylvest')
%
%       MaxNumComp: Maximum number of components for the vectors containig
%       the performance measure.
%       
%       Verbose: When set to 1, the progress of the computation is printed
%       (default:0).


OptionNames=['Mode      ';
             'MaxNumComp';
             'Verbose   '];
OptionTypes=['char   ';
             'numeric';
             'numeric'];

OptionValues{1}=['Direct'; 
                 'Sylves'];
 
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.Mode='Sylves';
options.MaxNumComp = 1000;
options.Verbose = 0;

% Parse Parameters
K=size(alpha,2);
Q_CT_MMAPK_ParsePara(D0,'D0',D,'D')
for i=1:K
    Q_CT_PH_ParsePara(alpha{i}, ['alpha_' int2str(i)], S{i}, ['S_', int2str(i)] )
end

% Parse Optional Parameters
options=Q_ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);


% Determine constants
m=size(D0,1);
smk=zeros(1,K+1);
smk(1)=0;
for i=1:K
    smk(i+1)=smk(i)+size(alpha{i},2);
end    
mser=smk(K+1);
mtot=mser*m;

% Test the load of the queue
Dsum=D0;
for i=1:K
    Dsum=Dsum+D{i};
end    
invD0=(-D0)^(-1);
pi=stat(eye(m)-Dsum/max(-diag(Dsum)));
lambdas=zeros(1,K);
mus=zeros(1,K);
lambda=sum(pi*(Dsum-D0));
for i=1:K
   lambdas(i)=sum(pi*D{i});
   temp=S{i}-sum(S{i},2)*alpha{i};
   if (size(S{i},1)>1)
        beta{i}=stat(eye(size(S{i},1))-temp/max(-diag(temp)));
   else
        beta{i}=1;
   end     
   mus(i)=-beta{i}*sum(S{i},2);
end   
load=lambdas*(1./mus)';
if load >= 1
    error('MATLAB:Q_CT_MMAPK_PHK_1:LoadExceedsOne',...
                        'The load %d of the system exceeds one',load);
end  

% Construct building blocks
LM=zeros(mtot,mtot);
Tser=zeros(mser,mser);
tser=zeros(mser,1);

for i=1:K
    Tser(smk(i)+1:smk(i+1),smk(i)+1:smk(i+1))=S{i};
    tser(smk(i)+1:smk(i+1),1)=-sum(S{i},2);
end
for i=1:K
    LM=LM+kron(D{i},tser*[zeros(1,smk(i)) alpha{i} zeros(1,smk(K+1)-smk(i+1))]);
end    

% Compute T iteratively
Told=zeros(mtot,mtot);
eyeTser=kron(eye(m),Tser);
D0eye=kron(D0,eye(mser));
Tnew=eyeTser;
if (options.Mode=='Direct')
    eyeD0eye=kron(eye(mtot),D0eye);
    while(max(max(abs(Told-Tnew)))>10^(-10))
        Told=Tnew;
        L=-reshape(eye(mtot),1,mtot^2)*(kron(Tnew',eye(mtot))+...
            eyeD0eye)^(-1);
        Tnew=eyeTser+reshape(L,mtot,mtot)'*LM;
    end
    L=reshape(L,mtot,mtot)';
else
    [U,Tr]=schur(D0,'complex');
    %U'*Tr*U = A
    U=kron(U,eye(mser));
    Tr=kron(Tr,eye(mser));
    while(max(max(abs(Told-Tnew)))>10^(-10))
        Told=Tnew;
        L=Q_Sylvest(U,Tr,Tnew);
        Tnew=eyeTser+L*LM;
    end
end


%Compute Witing Time and Sojourn Time distributions

% Compute S
theta_tot=zeros(1,mtot);
for i=1:K
    theta_tot=theta_tot+kron(pi*D{i},[zeros(1,smk(i)) beta{i} zeros(1,smk(K+1)-smk(i+1))]/mus(i));
end
theta_tot=theta_tot/load;
nonz=find(theta_tot>0);
theta_tot_red=theta_tot(nonz);
Tnewr=Tnew(:,nonz);
Tnewr=Tnewr(nonz,:);
Smat=diag(theta_tot_red)^(-1)*Tnewr'*diag(theta_tot_red);

% alpha vector of PH representation of Sojourn time
for i=1:K
    temp=zeros(mser,1);
    temp(smk(i)+1:smk(i+1))=tser(smk(i)+1:smk(i+1));
    soj_alpha{i}=load*(diag(theta_tot)*kron(ones(m,1),temp))'/lambdas(i);
    soj_alpha{i}=soj_alpha{i}(nonz);
end
soj_alpha{K+1}=load*(diag(theta_tot)*kron(ones(m,1),tser))'/lambda;
soj_alpha{K+1}=soj_alpha{K+1}(nonz);

% alpha vector of PH representation of Waiting time
for i=1:K
    wait_alpha{i}=load*(diag(theta_tot)*L*kron(sum(D{i},2),tser))'/lambdas(i);
    wait_alpha{i}=wait_alpha{i}(nonz);
end
wait_alpha{K+1}=load*(diag(theta_tot)*L*kron(sum(Dsum-D0,2),tser))'/lambda;
wait_alpha{K+1}=wait_alpha{K+1}(nonz);


% Overall queue length distribution
[V,NBAR]=hess(Tnew);
[U,Tr]=schur(D0,'complex');
%U'*Tr*U = D0
U=kron(U,eye(mser));
Tr=kron(Tr,eye(mser));

n=1;
Cn=Tnew;
qlT=1-load;
while (sum(qlT)<1-10^(-10) && n < 1+options.MaxNumComp)
    F=-V'*Cn*U;
    tempmat=zeros(mtot,mtot-1);
    for k=1:mtot
        if (k==1)
            temp=F(:,k);
        else
            temp=F(:,k)-Ln(:,1:k-1)*Tr(1:k-1,k);
        end
        Ln(:,k)=(NBAR+eye(mtot)*Tr(k,k))\temp;
    end
    Ln=real(V*Ln*U');
    Cn=Ln*kron(Dsum-D0,eye(mser));
    qlT(n+1)=-load*theta_tot*sum(Ln,2);
    n=n+1;
end
if (n == 1+options.MaxNumComp)
    warning('Maximum Number of Components %d reached',n-1);
end


% Queue length distribution per type
% Remark:
% L(0) solution to TL(0)+L(0)(Dsum-Dk)=-I
% L(n) to TL(n)+L(n)(Dsum-Dk)=-L(n-1)*kron(Dk,eye(mser))
[V,NBAR]=hess(Tnew);
for i=1:K
    ek=zeros(mser,1);
    ek(smk(i)+1:smk(i+1),1)=ones(smk(i+1)-smk(i),1);
    ek=kron(ones(m,1),ek);
    nek=ones(mtot,1)-ek;
    [U,Tr]=schur(Dsum-D{i},'complex');
    %U'*Tr*U = Dsum-D{i}
    U=kron(U,eye(mser));
    Tr=kron(Tr,eye(mser));
    
    ql{i}(1)=1-load; 
    n=1;
    Cn=Tnew;
    while (sum(ql{i})<1-10^(-10) && n < 1+options.MaxNumComp)
        F=-V'*Cn*U;
        tempmat=zeros(mtot,mtot-1);
        for k=1:mtot
            if (k==1)
                temp=F(:,k);
            else
                temp=F(:,k)-Ln(:,1:k-1)*Tr(1:k-1,k);
            end
            Ln(:,k)=(NBAR+eye(mtot)*Tr(k,k))\temp;
        end
        Ln=real(V*Ln*U');
        Cn=Ln*kron(D{i},eye(mser));
        ql{i}(n)=ql{i}(n)-load*theta_tot*Ln*nek;
        ql{i}(n+1)=-load*theta_tot*Ln*ek;
        n=n+1;
    end
    if (n == 1+options.MaxNumComp)
        warning('Maximum Number of Components %d reached',n-1);
    end
end