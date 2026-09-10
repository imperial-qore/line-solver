function [G,R,U]=MG1_IS(D,varargin)
%MG1_IS Invariant Subspace for M/G/1-Type Markov Chains [Akar,Sohraby] 
%
%   G=MG1_IS(A) computes the minimal nonnegative solution to the 
%   matrix equation G = A0 + A1 G + A2 G^2 + A3 G^3 + ... + A_max G^max, 
%   where A = [A0 A1 A2 A3 ... A_max] has m rows and m*max columns and is
%   a nonnegative matrix, with (A0+A1+A2+...+A_max) irreducible and 
%   stochastic
%
%   Optional Parameters:
%   
%       MaxNumIt: Maximum number of iterations (default: 50)
%       Verbose: The residual error is printed at each step when set to 1,  
%                (default:0)
%       Mode: 'MSignStandard' uses the matrix sign approach 
%             'MSignBalzer' uses the matrix sign approach with Balzer
%               acceleration
%             'Schur' relies on a Schur decomposition       
%             (default: MSignBalzer)

OptionNames=['Mode       ';
             'Verbose    ';
             'MaxNumIt   '];
OptionTypes=['char   ';
             'numeric';
             'numeric'];

OptionValues{1}=['MSignStandard';
                 'MSignBalzer  ';
                 'Schur        '];         
             
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.Mode='MSignBalzer';
options.Verbose=0;
options.MaxNumIt=50;

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

% check whether G is known explicitly
G=MG1_EG(D,options.Verbose);
if (~isempty(G))
    return
end

epsilon=10^(-14);
m=size(D,1);
f=size(D,2)/m-1;

maxd=size(D,2)/m-1;
sumD=D(:,maxd*m+1:end);
beta=sum(sumD,2);
% beta = (D_maxd)e + (D_maxd + D_maxd-1)e + ... + (Dmaxd+...+D1)e
for i=maxd-1:-1:1
    sumD=sumD+D(:,i*m+1:(i+1)*m);
    beta=beta+sum(sumD,2);
end
sumD=sumD+D(:,1:m);
theta=stat(sumD);
drift=theta*beta-1;

% Step 1, F(z)=z-A(z)
for i=0:f
    F{i+1}=-D(:,i*m+1:(i+1)*m);
end
F{2}=eye(m)+F{2};

% Step 2, H(s)=sum_{i=0}^f F_i (1-s)^(f-i) (1+s)^i = sum_{i=0}^f H_i s^i
for i=0:f
    H{i+1}=zeros(m,m);
end
for i=0:f
    temp1=[1 -1];
    temp2=[1 1];
    con1=1;
    con2=1;
    for j=1:f-i
        con1=conv(con1,temp1);
    end
    for j=1:i
        con2=conv(con2,temp2);
    end
    contrib=conv(con1,con2);
    for j=0:f
        H{j+1}=H{j+1}+contrib(j+1)*F{i+1};
    end    
end    

clear F;

% Step 3, \hat{H}_i = H_f^-1*H_i
H{f+1}=inv(H{f+1});
for i=0:f-1
    hatH{i+1}=H{f+1}*H{i+1};
end

% Step 4, y, xT
y=[ones(m,1); zeros(m*(f-1),1)];
x0T=[zeros(1,m) 1] / [hatH{1} ones(m,1)];
for i=1:f-1
    xT(1,(i-1)*m+1:i*m)=x0T*hatH{i+1};
end
xT(1,(f-1)*m+1:f*m)=x0T;

% Step 5, E_m in Zold
Zold=zeros(m*f,m*f);
for i=1:f-1
    Zold((i-1)*m+1:i*m,i*m+1:(i+1)*m)=eye(m);
end
for i=0:f-1
    Zold(m*(f-1)+1:m*f,i*m+1:(i+1)*m)=-hatH{i+1};
end
y=y/(xT*y);
Zold=Zold+sign(drift)*y*xT;

% Step 6, classic matrix sign function algorithm
if ( exist('ordschur') ~= 5 || ~strcmp(options.Mode,'Schur'))
    % Step 6, classic matrix sign function algorithm
    if (strcmp('Schur',options.Mode))
        fprintf('Ordschur not supported by current MATLAB version\n');
        fprintf('An automatic switch is performed to the MSignBalzer Mode\n');
        drawnow;
    end
    Znew=(Zold+inv(Zold))/2;
    numit=0;
    check=1;
    while (check > epsilon && numit < options.MaxNumIt)
        numit=numit+1;
        Zold=Znew;
        if (strcmp(options.Mode,'MSignStandard'))
            determ=1/2;
        else
            determ=(1+abs(det(Zold))^(1/(m*f)))^(-1);
        end
        Znew=determ*Zold+(1-determ)*inv(Zold);
        check=norm(Znew-Zold,1)/norm(Zold,1);
        if (options.Verbose==1)
            fprintf('Check after %d iterations: %d\n',numit,check);
            drawnow;
        end
    end
    if (numit == options.MaxNumIt && check > epsilon)
        warning('Maximum Number of Iterations %d reached: T may not have m columns',numit);
    end
    % Step 7,
    T=orth(Znew-eye(m*f));
else
    [T,U]=schur(Zold);
    [T,U]=ordschur(T,U,'lhp');
    T=T(:,1:m);
end    

% Step 8,
G=(T(1:m,:)+T(m+1:2*m,:))*inv(T(1:m,:)-T(m+1:2*m,:));
if (options.Verbose==1)
    temp=D(:,end-m+1:end);
    for i=size(D,2)/m-1:-1:1
        temp=D(:,(i-1)*m+1:i*m)+temp*G;
    end
    res_norm=norm(G-temp,inf);
    fprintf('Final Residual Error for G: %d\n',res_norm);
end

