function G=NSF_GHT(A,N,varargin)
%NSF_GHT Gail, Hantler, Taylor algorithm for Non-Skip-Free Type Markov Chains 
%
%   G=NSF_GHT(A,N) computes the minimal nonnegative solution to the 
%   matrix equation G = C0 + C1 G + C2 G^2 + C3 G^3 + ... + C_max G^max,
%   of the reblocked Non-Skip-Free Markov chain.
%   A = [A0 A1 A2 A3 ... A_max] has m rows and m*(max+1) columns and is
%   a nonnegative matrix, with (A0+A1+A2+...+A_max) irreducible and 
%   stochastic. The parameter N identifies the number of non-zero blocks
%   below the main diagonal, the matrices Ci are of dimension mN. 
%
%   Optional Parameters:
%   
%       MaxNumIt: Maximum number of iterations (default: 10000)
%       Verbose: When set to k, the residual error is printed every 
%                k steps (default:0)
%       FirstBlockRow: When set to one, only the first block row of G
%                      is returned, which fully characterizes G (default:0)

OptionNames=['MaxNumIt     ';
             'Verbose      ';
             'FirstBlockRow'];
OptionTypes=['numeric';
             'numeric';
             'numeric'];
OptionValues=[];
 
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.MaxNumIt=10000;
options.Verbose=0;
options.FirstBlockRow=0;

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

m=size(A,1);
K=size(A,2)/m-1; % degree of A(z) 
numit=0;
check=1;
G=A(:,1:m*N);

while(check > 10^(-14) && numit < options.MaxNumIt)
    numit=numit+1;
    Gold=G;
    temp=G;
    G=A(:,1:m*N)+A(:,m*N+1:m*(N+1))*G;
    for j=N+1:K
        temp=[zeros(m) temp(:,1:(N-1)*m)]+temp(:,(N-1)*m+1:end)*Gold;
        G=G+A(:,j*m+1:(j+1)*m)*temp;
    end    
    check=norm(G-Gold,inf);
    if (~mod(numit,options.Verbose))
        fprintf('Check after %d iterations: %d\n',numit,check);
        drawnow;
    end
end

if (numit == options.MaxNumIt && check > 10^(-14))
    warning('Maximum Number of Iterations %d reached',numit);
end

if ( ~options.FirstBlockRow )
    % compute remaining blockrows of G
    G(m+1:N*m,1:N*m)=zeros(m*(N-1),m*N);
    for j=2:N
        G((j-1)*m+1:j*m,:)=[zeros(m) G((j-2)*m+1:(j-1)*m,1:(N-1)*m)]+...
            G((j-2)*m+1:(j-1)*m,(N-1)*m+1:end)*G(1:m,:);
    end    
end    

if (options.Verbose>0)
    % we compute the residual error of first blockrow
    if (~options.FirstBlockRow)
        A=[A zeros(m,(N-1-mod(K,N))*m)];
        Nb=size(A,2)/(m*N)-1;
        Gcheck=A(:,Nb*N*m+1:end);
        for j=Nb-1:-1:0
            Gcheck=A(:,j*N*m+1:(j+1)*N*m)+Gcheck*G;
        end
        res_norm=norm(G(1:m,:)-Gcheck,inf);
        fprintf('Final Residual Error for G: %d\n',res_norm);
    end
end
