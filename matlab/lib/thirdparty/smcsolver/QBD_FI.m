function [G,R,U]=QBD_FI(A0,A1,A2,varargin)
%QBD_FI Functional Iterations for Quasi-Birth-Death Markov Chains [Neuts]  
%
%   DISCRETE TIME CASE:
%
%   G=QBD_FI(A0,A1,A2) computes the minimal nonnegative solution to the 
%   matrix equation G = A0 + A1 G + A2 G^2, where A,B and C are square 
%   nonnegative matrices, with (A0+A1+A2) irreducible and stochastic
%
%   [G,R]=QBD_FI(A0,A1,A2) also provides the minimal nonnegative solution 
%   to the matrix equation R = A2 + R A1 + R^2 A0
%
%   [G,R,U]=QBD_FI(A0,A1,A2) also provides the minimal nonnegative solution 
%   to the matrix equation U = A1 + A2 (I-U)^(-1) A0
%
%   CONTINUOUS TIME CASE:
%
%   G=QBD_FI(A0,A1,A2) computes the minimal nonnegative solution to the 
%   matrix equation 0 = A0 + A1 G + A2 G^2, where A,B and C are square 
%   nonnegative matrices, with (A0+A1+A2) having row sums equal to zero 
%
%   [G,R]=QBD_FI(A0,A1,A2) also provides the minimal nonnegative solution 
%   to the matrix equation 0 = A2 + R A1 + R^2 A0
%
%   [G,R,U]=QBD_FI(A0,A1,A2) also provides the minimal nonnegative solution 
%   to the matrix equation U = A1 + A2 (-U)^(-1) A0
%
%   Optional Parameters:
%   
%       MaxNumIt: Maximum number of iterations (default: 10000)
%       Mode: 'Traditional': G(n+1) = (I-A1)^(-1) * (A0 + A2 * G^2)
%             'Natural': G(n+1) = A0 + (A1 + A2*G(n))*G(n)
%             'U-Based': G(n+1) = (I-A1-A2*G(n))^(-1)*A0
%             'Shift<Mode>': where <Mode> is Traditional, Natural or
%             U-Based uses the Shift Technique
%             (default:'U-based')    
%       Verbose: When set to k, the residual error is printed every 
%                k steps (default:0)
%       StartValue: Starting value for iteration (default: 0)

OptionNames=['Mode       '; 
             'MaxNumIt   ';
             'Verbose    ';
             'StartValue '];
OptionTypes=['char   '; 
             'numeric';
             'numeric';
             'numeric'];
OptionValues{1}=['Traditional      '; 
                 'Natural          ';
                 'U-Based          ';
                 'ShiftTraditional ';
                 'ShiftNatural     ';
                 'ShiftU-Based     ';];
 
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.Mode='U-Based';
options.MaxNumIt=10000;
options.Verbose=0;
m=size(A1,1);
options.StartValue=zeros(m,m);

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

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

% check whether G is known explicitly
[G,R,U]=QBD_EG(A0,A1,A2,options.Verbose,nargout);
if (~isempty(G))
    return
end

numit=0;
check=1;
G=options.StartValue;

% Shift Technique
if (strfind(options.Mode,'Shift')>0)
    theta=statvec(A0+A1+A2);
    drift=theta*sum(A0,2)-theta*sum(A2,2);
    if (drift < 0) % MC is transient -> use the dual MC
        if (nargout > 1 | options.Verbose>0)
            A2old=A2;
        end
        A2=A2-ones(m,1)*(theta*A2);
        A1=A1+ones(m,1)*(theta*A0);
    else
        uT=ones(1,m)/m;
        A1old=A1;
        if (nargout > 2 | options.Verbose>0) % store A0old to compute U
            A0old=A0;
        end
        A0=A0-sum(A0,2)*uT;
        A1=A1+sum(A2,2)*uT;
    end
end    

if (strfind(options.Mode,'Natural')>0)
    while(check > 10^(-14) & numit < options.MaxNumIt)
        Gold=G;
        G=(A2*G+A1)*G+A0;
        check=norm(G-Gold,inf);
        numit=numit+1;
        if (~mod(numit,options.Verbose))
            fprintf('Check after %d iterations: %d\n',numit,check);
            drawnow;
        end
    end   
end    

if (strfind(options.Mode,'Traditional')>0)
    invA1=(eye(m)-A1)^(-1);
    while(check > 10^(-14) & numit < options.MaxNumIt)
        Gold=G;
        G=invA1*(A0+A2*G^2);
        check=norm(G-Gold,inf);
        numit=numit+1;
        if (~mod(numit,options.Verbose))
            fprintf('Check after %d iterations: %d\n',numit,check);
            drawnow;
        end
    end   
end    

if (strfind(options.Mode,'U-Based')>0)
    while(check > 10^(-14) & numit < options.MaxNumIt)
        Gold=G;
        G=(eye(m)-A1-A2*G)^(-1)*A0;
        check=norm(G-Gold,inf);
        numit=numit+1;
        if (~mod(numit,options.Verbose))
            fprintf('Check after %d iterations: %d\n',numit,check);
            drawnow;
        end
    end   
end  
if (numit == options.MaxNumIt)
    warning('Maximum Number of Iterations %d reached',numit);
end

if (strfind(options.Mode,'Shift')>0)
    if (drift < 0) % transient
        if (nargout > 1 | options.Verbose >0)
            A1=A1-ones(m,1)*theta*A0; % restore original A1
            A2=A2old; % restore original A2
        end     
    else % pos recurrent
        G=G+ones(m,1)*uT;
        if (nargout > 1 | options.Verbose >0)
            A1=A1-sum(A2,2)*uT; % restore original A1 
        end
        if (nargout > 2 | options.Verbose >0)
            A0=A0old; % restore original A0
        end
    end
end    

if (options.Verbose>0)
    res_norm=norm(G-A0-(A1+A2*G)*G,inf);
    fprintf('Final Residual Error for G: %d\n',res_norm);
end

% Compute R
if (nargout > 1)
    R=A2*(eye(m)-(A1+A2*G))^(-1);
    if (options.Verbose>0)
        res_norm=norm(R-A2-R*(A1+R*A0),inf);
        fprintf('Final Residual Error for R: %d\n',res_norm);
    end
end    

% Compute U
if (nargout > 2)
    U=A1+R*A0;
    if (options.Verbose>0)
        res_norm=norm(U-A1-A2*(eye(m)-U)^(-1)*A0,inf);
        fprintf('Final Residual Error for U: %d\n',res_norm);
    end
    if (continues)
        U=lamb*(U-eye(m));
    end    
end    