function [G,R,U]=QBD_LR(A0,A1,A2,varargin)
%QBD_LR Logaritmic reduction for Quasi-Birth-Death Markov Chains [Latouche,Ramaswami] 
%
%   DISCRETE TIME CASE:
%
%   G=QBD_LR(A0,A1,A2) computes the minimal nonnegative solution to the 
%   matrix equation G = A0 + A1 G + A2 G^2, where A,B and C are square 
%   nonnegative matrices, with (A0+A1+A2) irreducible and stochastic
%
%   [G,R]=QBD_LR(A0,A1,A2) also provides the minimal nonnegative solution 
%   to the matrix equation R = A2 + R A1 + R^2 A0
%
%   [G,R,U]=QBD_LR(A0,A1,A2) also provides the minimal nonnegative solution 
%   to the matrix equation U = A1 + A2 (I-U)^(-1) A0
%
%   CONTINUOUS TIME CASE:
%
%   G=QBD_LR(A0,A1,A2) computes the minimal nonnegative solution to the 
%   matrix equation 0 = A0 + A1 G + A2 G^2, where A,B and C are square 
%   nonnegative matrices, with (A0+A1+A2) having row sums equal to zero 
%
%   [G,R]=QBD_LR(A0,A1,A2) also provides the minimal nonnegative solution 
%   to the matrix equation 0 = A2 + R A1 + R^2 A0
%
%   [G,R,U]=QBD_LR(A0,A1,A2) also provides the minimal nonnegative solution 
%   to the matrix equation U = A1 + A2 (-U)^(-1) A0
%
%   Optional Parameters:
%   
%       MaxNumIt: Maximum number of iterations (default: 50)
%       Verbose: The residual error is printed at each step when set to 1,  
%                (default:0)
%       Mode: 'Basic' uses the Basic Cyclic Reduction 
%             'Shift' uses the shift technique to accelarate convergence
%             (default: 'Shift')  
%       RAPComp: set to 1 if the QBD has RAP components

OptionNames=['Mode       '; 
             'MaxNumIt   ';
             'Verbose    ';
             'RAPComp    '];
OptionTypes=['char   '; 
             'numeric';
             'numeric';
             'numeric'];
OptionValues{1}=['Basic'; 
                 'Shift'];
 
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.Mode='Shift';
options.MaxNumIt=50;
options.Verbose=0;
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

% Start LR
m=size(A1,1);


% Shift technique
if (options.Mode=='Shift')
    theta=stat(A0+A1+A2);
    drift=theta*sum(A0,2)-theta*sum(A2,2);
    if (drift < 0) % MC is transient -> use the dual MC
        if (nargout > 1 | options.Verbose==1)
            A2old=A2;
        end
        A2=A2-ones(m,1)*(theta*A2);
        A1=A1+ones(m,1)*(theta*A0);
    else
        uT=ones(1,m)/m;
        if (nargout > 2 | options.Verbose==1) % store A0old to compute U
            A0old=A0;
        end
        A0=A0-sum(A0,2)*uT;
        A1=A1+sum(A2,2)*uT;
    end
end    

% Start of Logaritmic Reduction (Basic)
B2=(eye(m) - A1)^(-1);
B0=B2*A2;
B2=B2*A0;
if (nargout <= 1 & options.Verbose ~= 1) % A1 and A2 only needed to compute R, U
    clear A1 A2;
end    
G=B2;
PI=B0;
check=1;
numit=0;
while(check > 10^(-14) & numit < options.MaxNumIt)
    A1star=B2*B0+B0*B2;
    A0star=B0^2;
    A2star=B2^2;
    B0=(eye(m)-A1star)^(-1);
    B2=B0*A2star;
    B0=B0*A0star;
    G=G+PI*B2;
    PI=PI*B0;
    check=min(norm(B0,inf),norm(B2,inf));
    numit=numit+1;
    if (options.Verbose==1)
        fprintf('Check after %d iterations: %d\n',numit,check);
        drawnow;
    end
end
if (numit == options.MaxNumIt)
    warning('Maximum Number of Iterations %d reached',numit);
end    
clear PI A0star A1star A2star B0 B2;

% Shift Technique
if (options.Mode=='Shift')
    if (drift < 0) % transient
        if (nargout > 1 | options.Verbose==1)
            A1=A1-ones(m,1)*theta*A0; % restore original A1
            A2=A2old; % restore original A2
        end     
    else % pos recurrent
        G=G+ones(m,1)*uT;
        if (nargout > 1 | options.Verbose==1)
            A1=A1-sum(A2,2)*uT; % restore original A1 
        end
        if (nargout > 2 | options.Verbose==1)
            A0=A0old; % restore original A0
        end
    end
end

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
