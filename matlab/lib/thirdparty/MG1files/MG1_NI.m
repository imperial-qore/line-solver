function G=MG1_NI(A,varargin)
% Newton Iteration for M/G/1-Type Markov Chains 
%
%   G=MG1_NI(A) solves G = A0 + A1 G + A2 G^2 + A3 G^3 + ... + A_max G^max, 
%   where A = [A0 A1 A2 A3 ... AN] has m rows and m*(N+1) columns and is
%   a nonnegative matrix, with (A0+A1+A2+...+AN irreducible and stochastic
%
%   Optional Parameters:
%
%       Mode: method used to solve the linear system at each iteration and
%       optional shifting
%               'DirectSum':        solves the system directly
%               'DirectSumShift':   solves the system directly + Shift
%               'RealSchur':        applies a real Schur decomposition (default)
%               'RealSchurShift':   applies a real Schur decomposition + Shift
%               'ComplexSchur':     applies a complex Schur decomposition
%               'ComplexSchurShift': applies a complex Schur decomposition + Shift
%       MaxNumIt: Maximum number of iterations (default: 50)
%       ShiftType: 'one' (default if shift mode is selected)
%                  'tau'
%                  'dbl'
%       EpsilonValue: Required accuracy, used as stop criterium (default:10^(-14)) 
%       Verbose: The residual error is printed at each step when set to 1,  
%                (default:0)


OptionNames=['Mode        '; 
             'MaxNumIt    ';
             'ShiftType   ';
             'EpsilonValue';
             'Verbose     '];
OptionTypes=['char   '; 
             'numeric';
             'char   ';
             'numeric';
             'numeric'];
OptionValues{1}=['DirectSum        ';
                 'DirectSumShift   ';
                 'RealSchur        ';
                 'RealSchurShift   ';
                 'ComplexSchur     ';
                 'ComplexSchurShift'];

OptionValues{3}=['one';
                 'tau';
                 'dbl'];

             
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.Mode='RealSchurShift';
options.MaxNumIt=50;
options.ShiftType='one';
options.Verbose=0;
options.EpsilonValue=10^(-14);

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

% check whether G is known explicitly
G=MG1_EG(A,options.Verbose);
if (~isempty(G))
    return
end

if (strfind(options.Mode,'Shift')>0)
    if (options.Verbose == 1)
        Aold=A;
    end
    [A,drift,tau,v]=MG1_Shifts(A,options.ShiftType);  
end


m = size(A,1);
N = size(A,2)/m-1;

G = zeros(m);
check = 1;
numit = 0;
while(check > options.EpsilonValue && numit < options.MaxNumIt)
    Gold = G;
     
    %Compute B: matrices premultiplying Y_k, and C: right-hand side
    B = zeros(m, (N+1)*m);
    B(:,N*m+1:(N+1)*m) = A(:,N*m+1:(N+1)*m);
    for i = N-1:-1:0
        B(:,i*m+1:(i+1)*m) = A(:,i*m+1:(i+1)*m) + B(:,(i+1)*m+1:(i+2)*m)*G;
    end
    C = G - B(:,1:m);
    B = B(:,m+1:end);
    B(:,1:m) = B(:,1:m) - eye(m);
    
    if strfind(options.Mode,'RealSchur') > 0
        Y = solveSylvPowersRealSchur_FW(G, B, C);
    elseif strfind(options.Mode,'ComplexSchur') > 0
        Y = solveSylvPowersComplexSchur_FW(G, B, C);
    else
        Y = solveSylvPowersDirectSum(G, B, C);
    end
    
    G = G + Y;
    check=norm(G-Gold,inf);
    numit = numit + 1;
    if (~mod(numit,options.Verbose))
        fprintf('Check after %d iterations: %d\n',numit,check);
    end
end
%numit
if (numit == options.MaxNumIt)
    warning('Maximum Number of Iterations %d reached',numit);
end  

if (strfind(options.Mode,'Shift'))
    switch options.ShiftType
        case 'one'
            G=G+(drift<1)*ones(m,m)/m;
        case 'tau'
            G=G+(drift>1)*tau*v*ones(1,m);
        case 'dbl'
            G=G+(drift<1)*ones(m,m)/m+(drift>1)*tau*v*ones(1,m);
    end    
end  

if (options.Verbose>0)
    if (strfind(options.Mode,'Shift'))
        A=Aold;
    end
    temp=A(:,end-m+1:end);
    for i=N:-1:1
        temp=A(:,(i-1)*m+1:i*m)+temp*G;
    end
    res_norm=norm(G-temp,inf);
    fprintf('Final Residual Error for G: %d\n',res_norm);
end