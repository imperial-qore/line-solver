function G=MG1_NI_LRAi(A0, Ahat, Gamma, varargin)
% Newton Iteration for M/G/1-Type Markov Chains, exploiting low-rank Ai,
% for i=1,...,N. 
%
%   G=MG1_NI_LRAi(A0, Ahat, Gamma) solves G = A0 + A1 G + A2 G^2 + A3 G^3 + ... + A_N G^N, 
%   where A0, A1,...,AN are mxm nonnegative matrices such that
%   (A0+A1+A2+...+AN) is irreducible and stochastic.
%   Ai can be written as Ai = Gamma*Aihat, for i = 1,...,N, with Aihat rxm
%   matrices and Gamma an mxr matrix.
%   Ahat = [A1hat A2hat A3hat ... ANhat] has r rows and m*N columns.
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
m = size(A0,1);
options.Mode='RealSchur';
options.MaxNumIt=50;
options.Verbose=0;
options.EpsilonValue=10^(-14);

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);


N = size(Ahat,2)/m;
r = size(Ahat,1);

check = 1;
numit = 1;

Uhat = zeros(r,m);
while(check > options.EpsilonValue && numit < options.MaxNumIt)
    Uhatold = Uhat;
    
    %Compute B: matrices premultiplying Y_k, and C: right-hand side
    Uinv=eye(m)+Gamma*inv(eye(r)-Uhatold*Gamma)*Uhatold;
    K = Uinv*A0;
    B = zeros(r, (N-1)*r);
    temp = Ahat(:,(N-1)*m+1:N*m);
    UinvG = Uinv*Gamma;
    for i=N-2:-1:0
        B(:,i*r+1:(i+1)*r) = temp*UinvG;
        temp = Ahat(:,i*m+1:(i+1)*m)+temp*K;
    end    
    C = Uhatold-temp;
    B = [-eye(r) B];
    
    checkB = 0;
    k = N;
    while k > 1 && checkB < options.EpsilonValue
        checkB = norm(B(:, (k-1)*r+1:k*r), 'inf');
        k = k-1;
    end
    
    if k < N-1
        B = B(:,1:(k+1)*r);
    end
    
    if strfind(options.Mode,'RealSchur') > 0
        Y = solveSylvPowersRealSchur_FW(K, B, C);
    elseif strfind(options.Mode,'ComplexSchur') > 0
        Y = solveSylvPowersComplexSchur_FW(K, B, C);
    else
        Y = solveSylvPowersDirectSum(K, B, C);
    end
    
    Uhat = Uhat + Y;
    check=norm(Uhat-Uhatold,inf);
    numit = numit + 1;
    if (~mod(numit,options.Verbose))
        fprintf('Check after %d iterations: %d\n',numit,check);
    end
end
G = inv(eye(m) - Gamma*Uhat)*A0;

if (numit == options.MaxNumIt)
    warning('Maximum Number of Iterations %d reached',numit);
end  

if (options.Verbose>0)
    A = A0;
    for j = 1:N
        A = [A Gamma*Ahat(:,(j-1)*m+1:j*m)];
    end
    temp=A(:,end-m+1:end);
    for i=N:-1:1
        temp=A(:,(i-1)*m+1:i*m)+temp*G;
    end
    res_norm=norm(G-temp,inf);
    fprintf('Final Residual Error for G: %d\n',res_norm);
end