function pi=MG1_pi(B,A,G,varargin)
%NSF_pi Stationary vector of a Non-Skip-Free Markov Chains [Neuts, Ramaswami] 
%
%   pi=NSF_pi(B,A,G) computes the stationary vector of a Non-Skip-Free
%   Markov chain with a transition matrix of the form
%
%               B10  B11  B12  B13  B14  ...
%               B20  B21  B22  B23  B24  ...
%               ...
%               BN0  BN1  BN2  BN3  BN4  ...
%               A0   A1   A2   A3   A4  ...
%       P  =    0    A0   A1   A2   A3  ...
%               0    0    A0   A1   A2 ...
%               ...
%
%   the input matrix G is the minimal nonnegative solution to the matrix 
%   equation G = C0 + C1 G + C2 G^2 + ... + Ccmax G^cmax of the
%   reblocked system, A = [A0 A1 ... Aamax] and B=[B00 B01 B02 ... B0bmax; 
%   B10 B11 B12 ... B1bmax; ... ; BN0 BN1 BN2 ... BNbmax]
%
%   pi=NSF_pi([],A,G) computes the stationary vector of an M/G/1-type
%   Markov chain with a transition matrix of the form, i.e., Bij = Aj
%   for i=1,...,N:
%
%               A0   A1   A2   A3   A4  ...          
%               A0   A1   A2   A3   A4  ...          
%               ...
%               A0   A1   A2   A3   A4  ...          
%               A0   A1   A2   A3   A4  ...
%       P  =    0    A0   A1   A2   A3  ...
%               0    0    A0   A1   A2 ...
%               ...
%
%   Optional Parameters:
%   
%       MaxNumComp: Maximum number of components (default: 500)
%       Verbose: The accumulated probability mass is printed at every 
%                n steps when set to n (default:0)
%       FirstBlockRow: When set to 1, it suffices to give the first
%                      blockrow of G as input (default:0)


OptionNames=['FirstBlockRow'
             'MaxNumComp   ';
             'Verbose      ';
             'Boundary     '];
OptionTypes=['numeric'
             'numeric';
             'numeric';
             'numeric'];
 
OptionValues=[];             
             
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.MaxNumComp=500;
options.Verbose=0;
options.Boundary=[];
options.FirstBlockRow=[];

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

m=size(A,1);
K=size(A,2)/m-1; % degree of A(z) 
N=size(G,2)/m;

if ( options.FirstBlockRow )
    % compute remaining blockrows of G
    G(m+1:N*m,1:N*m)=zeros(m*(N-1),m*N);
    for j=2:N
        G((j-1)*m+1:j*m,:)=[zeros(m) G((j-2)*m+1:(j-1)*m,1:(N-1)*m)]+...
            G((j-2)*m+1:(j-1)*m,(N-1)*m+1:end)*G(1:m,:);
    end    
end    

% construct reblocked Ci matrices
extraBlocksC=N-1-mod(K-1,N);
C=zeros(N*m,m*(N-1)+size(A,2)+m*extraBlocksC);
C(1:m,1:(K+1)*m)=A;
for i=1:N-1
    C(i*m+1:(i+1)*m,i*m+1:(i+K+1)*m)=A;
end

if (~isempty(B))
    extraBlocksB=N-1-mod(K,N);
    B=[B zeros(N*m,m*extraBlocksB)];
    pi=MG1_pi(B,C,G,'Verbose',options.Verbose,'MaxNumComp',options.MaxNumComp/N);
else
    FirstRowSumsG=G(1:m,1:m);
    for i=2:N
        FirstRowSumsG(1:m,(i-1)*m+1:i*m)=...
            FirstRowSumsG(1:m,(i-2)*m+1:(i-1)*m)+G(1:m,(i-1)*m+1:i*m);
    end
    ghat=stat(FirstRowSumsG(:,(N-1)*m+1:end));
    pi0=ghat*G(1:m,:);
    % to nomalize pi0 we need to compute the first m components of mu1
    g=ghat*FirstRowSumsG;
    beta=zeros(m*N,1);
    % beta = Cmax e + (Cmax + Cmax-1) e + ... (Cmax + ... + C1)e
    CSum=zeros(m*N,m*N);
    for i=size(C,2)/(m*N):-1:2
        CSum=CSum+C(:,(i-1)*m*N+1:i*m*N);
        beta=beta+sum(CSum,2);
    end    
    CSum=CSum+C(:,1:m*N);
    temp=(eye(m*N)-CSum+(ones(m*N,1)-beta)*g)^(-1)*ones(N*m,1);
    temp=([eye(m) zeros(m,(N-1)*m)]-G(1:m,:)+ones(m,1)*g)*temp;
    pi0=pi0/(ghat*temp);
    % compute B matrices of M/G/1
    B=kron(ones(N,1),A);
    extraBlocksB=N-1-mod(K,N);
    B=[B zeros(N*m,m*extraBlocksB)];
    pi=MG1_pi(B,C,G,'InputPiZero',pi0,'Verbose',options.Verbose,...
        'MaxNumComp',options.MaxNumComp/N);
end    