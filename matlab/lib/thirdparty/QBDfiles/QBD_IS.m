function [G,R,U]=QBD_IS(A0,A1,A2,varargin)
%QBD_IS Invariant Subspace for Quasi-Birth-Death Markov Chains [Akar, Sohraby] 
%
%   DISCRETE TIME CASE:
%
%   G=QBD_IS(A0,A1,A2) computes the minimal nonnegative solution to the 
%   matrix equation G = A0 + A1 G + A2 G^2, where A,B and C are square 
%   nonnegative matrices, with (A0+A1+A2) irreducible and stochastic
%
%   [G,R]=QBD_IS(A0,A1,A2) also provides the minimal nonnegative solution 
%   to the matrix equation R = A2 + R A1 + R^2 A0
%
%   [G,R,U]=QBD_IS(A0,A1,A2) also provides the minimal nonnegative solution 
%   to the matrix equation U = A1 + A2 (I-U)^(-1) A0
%
%   CONTINUOUS TIME CASE:
%
%   G=QBD_IS(A0,A1,A2) computes the minimal nonnegative solution to the 
%   matrix equation 0 = A0 + A1 G + A2 G^2, where A,B and C are square 
%   nonnegative matrices, with (A0+A1+A2) having row sums equal to zero 
%
%   [G,R]=QBD_IS(A0,A1,A2) also provides the minimal nonnegative solution 
%   to the matrix equation 0 = A2 + R A1 + R^2 A0
%
%   [G,R,U]=QBD_IS(A0,A1,A2) also provides the minimal nonnegative solution 
%   to the matrix equation U = A1 + A2 (-U)^(-1) A0
%
%   Optional Parameters:
%   
%       MaxNumIt: Maximum number of iterations (default: 50)
%       Verbose: The residual error is printed at each step when set to 1,  
%                (default:0)
%       Mode: 'MSignStandard' uses the matrix sign approach 
%             'MSignBalzer' uses the matrix sign approach with Balzer acceleration
%             'Schur' relies on an ordered Schur decomposition to find the 
%             invariant subspace (default: Schur)
%       RAPComp: set to 1 if the QBD has RAP components

OptionNames=['Mode       ';
             'Verbose    ';
             'MaxNumIt   ';
             'RAPComp    '];
OptionTypes=['char   ';
             'numeric';
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
options.Mode='Schur';
options.Verbose=0;
options.MaxNumIt=50;
options.RAPComp=0;

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

if(~options.RAPComp)
    % Convert to discrete time problem, if needed
    m=size(A1,1);
    continues=0;
    if (sum(diag(A1)<0)) % continues time
        continues=1;
        lamb=max(-diag(A1));
        A0=A0/lamb;
        A1=A1/lamb+eye(m);
        A2=A2/lamb;
    end

    % Parse Parameters
    QBD_ParsePara(A0,A1,A2);
else
    % Parse Parameters
    QBD_RAP_ParsePara(A0,A1,A2);
    
    % Convert to discrete time problem - uniformization
    m=size(A1,1);
    continues=1;
    lamb=max(-diag(A1));
    A0=A0/lamb;
    A1=A1/lamb+eye(m);
    A2=A2/lamb;
end

% check whether G is known explicitly
[G,R,U]=QBD_EG(A0,A1,A2,options.Verbose,nargout);
if (~isempty(G))
    return
end

epsilon=10^(-12);
f=2;
m=size(A1,1);


% Step 1, F(z)=zD(z)-N(z), with A(z) = D^(-1)(z)N(z)
% For QBD: A(z) is a polynomial matrix in z with degree f. Thus, D(z)=I.
theta=stat(A0+A1+A2);
drift=theta*sum(A0,2)-theta*sum(A2,2);
F{1}=-A0;
F{2}=eye(m)-A1;
F{3}=-A2;

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
Zold=Zold-sign(drift)*y*xT;

if ( exist('ordschur') ~= 5 || ~strcmp(options.Mode,'Schur'))
    % Step 6, classic matrix sign function algorithm
    if (strcmp('Schur',options.Mode))
        fprintf('Ordschur not supported by current MATLAB version\n');
        fprintf('An automatic switch is performed to the MSignBalzer Mode\n');
        drawnow;
    end
    numit=0;
    check=1;
    while (check > epsilon && numit < options.MaxNumIt)
        numit=numit+1;
        if (strcmp(options.Mode,'MSignStandard'))
            determ=1/2;
        else % options.Mode = 'MSignBalzer' or switched from 'Schur' Mode
            determ=(1+abs(det(Zold))^(1/(m*f)))^(-1);
            determ=min(determ,1-10^(-3));
        end    
        Znew=determ*Zold+(1-determ)*inv(Zold);
        check=norm(Znew-Zold,1)/norm(Zold,1);
        if (options.Verbose==1)
            fprintf('Check after %d iterations: %d\n',numit,check);
            drawnow;
        end
        Zold=Znew;
    end
    if (numit == options.MaxNumIt && check > epsilon)
        warning('Maximum Number of Iterations %d reached: T may not have m columns',numit);
    end
    % Step 7, 
    T=orth(Znew-eye(m*f));
else    
    % Theorem 7 (Akar, Oguz, Sohraby -- SM 2000)
    [T,D]=schur(Zold);
    [T,~]=ordschur(T,D,'lhp');
    T=T(:,1:m);
end

% Step 8,
G=(T(1:m,:)+T(m+1:2*m,:))*inv(T(1:m,:)-T(m+1:2*m,:));

if (options.Verbose==1)
    res_norm=norm(G-A0-(A1+A2*G)*G,inf);
    fprintf('Final Residual Error for G: %d\n',res_norm);
end

% Compute R
if (nargout > 1)
    R=A2*(eye(m)-(A1+A2*G))^(-1);
    if (options.Verbose==1)
        res_norm=norm(R-A2-R*(A1+R*A0),inf);
        fprintf('Final Residual Error for R: %d\n',res_norm);
    end
end   

% Compute U
if (nargout > 2)
    U=A1+R*A0;
    if (options.Verbose==1)
        res_norm=norm(U-A1-A2*(eye(m)-U)^(-1)*A0,inf);
        fprintf('Final Residual Error for U: %d\n',res_norm);
    end
    if (continues)
        U=lamb*(U-eye(m));
    end    
end    

