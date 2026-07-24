function [G,R,U]=QBD_NI(A0,A1,A2,varargin)
%QBD_NI Newton Iteration for Quasi-Birth-Death Markov Chains [Ramaswami, 
%   Latouche,Bini,Meini] 
%
%   DISCRETE TIME CASE:
%
%   G=QBD_NI(A0,A1,A2) computes the minimal nonnegative solution to the 
%   matrix equation G = A0 + A1 G + A2 G^2, where A,B and C are square 
%   nonnegative matrices, with (A0+A1+A2) irreducible and stochastic
%
%   [G,R]=QBD_NI(A0,A1,A2) also provides the minimal nonnegative solution 
%   to the matrix equation R = A2 + R A1 + R^2 A0
%
%   [G,R,U]=QBD_NI(A0,A1,A2) also provides the minimal nonnegative solution 
%   to the matrix equation U = A1 + A2 (I-U)^(-1) A0
%
%   CONTINUOUS TIME CASE:
%
%   G=QBD_NI(A0,A1,A2) computes the minimal nonnegative solution to the 
%   matrix equation 0 = A0 + A1 G + A2 G^2, where A,B and C are square 
%   nonnegative matrices, with (A0+A1+A2) having row sums equal to zero 
%
%   [G,R]=QBD_NI(A0,A1,A2) also provides the minimal nonnegative solution 
%   to the matrix equation 0 = A2 + R A1 + R^2 A0
%
%   [G,R,U]=QBD_NI(A0,A1,A2) also provides the minimal nonnegative solution 
%   to the matrix equation U = A1 + A2 (-U)^(-1) A0
%
%   Optional Parameters:
%   
%       MaxNumIt: Maximum number of iterations (default: 50)
%       Verbose: The residual error is printed at each step when set to 1,  
%                (default:0)
%       Mode: 'Sylvest' solves a Sylvester matrix equation at each step 
%             using an Hessenberg algorithm
%             'Estimat' estimates the solution of the Sylvester equation
%             'DirectSum' solves the Sylvester matrix equation at each
%             step by rewriting it as a (large) system of linear equations
%             (default: 'Sylvest')
%       RAPComp: set to 1 if the QBD has RAP components

OptionNames=[
%             'ProgressBar'; 
             'Mode       '; 
             'MaxNumIt   ';
             'Verbose    ';
             'RAPComp    '];
OptionTypes=[
%             'numeric'; 
             'char   '; 
             'numeric';
             'numeric';
             'numeric'];
OptionValues{1}=['Sylvest  '; 
                 'Estimat  ';
                 'DirectSum'];
 
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
%options.ProgressBar=0;
options.Mode='Sylvest';
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

% Start NI
m=size(A1,1);

R=zeros(m,m);
check=1;
numit=0;
while (check > 10^(-12) && numit < options.MaxNumIt)
    numit=numit+1;
    if (numit==1)%R=0;
        Yk=A2*(eye(m)-A1)^(-1);
    else
        if (strcmp(options.Mode,'Estimat'))
            FRk=(A2+R*(A1-eye(m)+R*A0)); % FRk = (A2+RA1+R^2A0)-R
            Zk=FRk*(eye(m)-A1)^(-1);
            Yk=FRk+Zk*A1+(R*Zk+Zk*R)*A0;
        else    
            D=-(A2+R*(A1-eye(m)+R*A0)); % D = R-(A2+RA1+R^2A0)
            C=A1+R*A0-eye(m);
            % Solve R*Yk*A0+Yk*C=D
            if (strcmp(options.Mode,'Sylvest'))
                Yk=QBD_NI_Sylvest(A0',R',C',D')';
            else
                Yk=(kron(A0',R)+kron(C',eye(m)))^(-1)*reshape(D,m^2,1);
                Yk=reshape(Yk,m,m);
            end    
        end    
    end
    R=R+Yk;
    check=norm(Yk,inf);
    if (options.Verbose==1)
        fprintf('Check after %d iterations: %d\n',numit,check);
        drawnow;
    end
end  
clear C D Yk;
if (numit == options.MaxNumIt && check > 10^(-12))
    warning('Maximum Number of Iterations %d reached',numit);
end    

% Compute G
G=(eye(m)-(A1+R*A0))^(-1)*A0;

if (options.Verbose==1)
    res_norm=norm(G-A0-(A1+A2*G)*G,inf);
    fprintf('Final Residual Error for G: %d\n',res_norm);
end

% R
if (nargout > 1)
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
