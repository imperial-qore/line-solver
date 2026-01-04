function G=MG1_NI_LRA0(A, A0hat, Gamma, varargin)
% Newton Iteration for M/G/1-Type Markov Chains, exploiting low-rank A0
%
%   G=MG1_NI_LRA0(A, A0hat, Gamma) solves  G = A0 + A1 G + ... + A_N G^N, 
%   where A = [A1 A2 A3 ... AN] has m rows and m*N columns and is
%   a nonnegative matrix. A0 is also nonnegative and can be written as 
%   A0 = A0hat*Gamma, where A0hat and Gamma are mxr and rxm matrices,
%   respectively. (A0+A1+A2+...+AN) is irreducible and stochastic.
%
%   Optional Parameters:
%
%       Mode: method used to solve the linear system at each iteration 
%               'DirectSum':        solves the system directly
%               'RealSchur':        applies a real Schur decomposition (default)
%               'ComplexSchur':     applies a complex Schur decomposition
%       MaxNumIt: Maximum number of iterations (default: 50)
%       EpsilonValue: Required accuracy, used as stop criterium (default:10^(-14)) 
%       Verbose: The residual error is printed at each step when set to 1,  
%                (default:0)


OptionNames=['Mode        '; 
             'MaxNumIt    ';
             'EpsilonValue';
             'Verbose     '];
OptionTypes=['char   '; 
             'numeric';
             'numeric';
             'numeric'];
OptionValues{1}=['DirectSum        ';
                 'RealSchur        ';
                 'ComplexSchur     '];

             
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.Mode='RealSchur';
options.MaxNumIt=50;
options.Verbose=0;
options.EpsilonValue=10^(-14);

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

% check whether G is known explicitly
G=MG1_EG([A0hat*Gamma A],options.Verbose);
if (~isempty(G))
    return
end

m = size(A,1);
N = size(A,2)/m;
r = size(Gamma,1);

check = 1;
numit = 1;

Ghat = zeros(m,r);
while(check > options.EpsilonValue && numit < options.MaxNumIt)
    Gold = Ghat;
   
    %Compute B: matrices premultiplying Y_k, and C: right-hand side
    B = zeros(m, (N+1)*r);
    B(:,(N-1)*r+1:N*r) = A(:,(N-1)*m+1:N*m)*Ghat;
    K = Gamma*Ghat;
    for i = N-2:-1:0
        B(:,i*r+1:(i+1)*r) = A(:,i*m+1:(i+1)*m)*Ghat + B(:,(i+1)*r+1:(i+2)*r)*K;
    end
    C = Ghat - B(:,1:r) - A0hat;
    B = reshape( [reshape(B(:,r+1:end),r*m,N); zeros((m-r)*m,N)], m,N*m);
    for i = 0:N-1
        B(:,i*m+1:(i+1)*m) = B(:,i*m+1:i*m+r)*Gamma;
    end
    %B = A(:,m+1:end) + B;
    B = A + B;
    B(:,1:m) = B(:,1:m) - eye(m);
   
    if strfind(options.Mode,'RealSchur') > 0
        Y = solveSylvPowersRealSchur_FW(K, B, C);
    elseif strfind(options.Mode,'ComplexSchur') > 0
        Y = solveSylvPowersComplexSchur_FW(K, B, C);
    else
        Y = solveSylvPowersDirectSum(K, B, C);
    end
    
    Ghat = Ghat + Y;
    check=norm(Ghat-Gold,inf);
    numit = numit + 1;
    if (~mod(numit,options.Verbose))
        fprintf('Check after %d iterations: %d\n',numit,check);
    end
end
G = Ghat*Gamma;

if (numit == options.MaxNumIt)
    warning('Maximum Number of Iterations %d reached',numit);
end  

if (options.Verbose>0)
    A = [A0hat*Gamma A];
    temp=A(:,end-m+1:end);
    for i=N:-1:1
        temp=A(:,(i-1)*m+1:i*m)+temp*G;
    end
    res_norm=norm(G-temp,inf);
    fprintf('Final Residual Error for G: %d\n',res_norm);
end