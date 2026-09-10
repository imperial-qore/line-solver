function Q_DT_SMK_ParsePara(C, C_name, processType)
% Q_DT_SMK_Service_ParsePara checks the validity of the input cell C as a
%   representation of a Discrete-Time marked Semi-Markov process.
%
%   If the SMK is a service process
%   C{k}=[Ck_1 Ck_2...Ck_tmax(k)], for k=1...K, where the entry (j,j') of
%   the mxm matrix Dk_i holds the probability that the service time of 
%   a type-k customer takes i time slots (i=1:tmax(k)), while the 
%   underlying phase changes from j to j'
%
%   If the SMK is an arrival process
%   C{k}=[Ck_1 Ck_2...Ck_tmax(k)], for k=1...K, where the entry (j,j') of
%   the mxm matrix Ck_i, holds the probabilities of having a type k
%   arrival, with an interarrival time of i time slots (i=1:tmax(k)), 
%   while the underlying phase changes from j to j'
%
%   processType: equals 1 if the SM[K] process is an arrival process and 0
%   if it represents a service process


if nargin < 3
    error('MATLAB:Q_DT_SMK_ParsePara:InvalidInput',...
            'The type of the SM[K] process (arrival or service) has to be specified');
end

% check numeric
K=size(C,2);
for i=1:K
    if (~isnumeric(C{i}))
        error('MATLAB:Q_DT_SMK_ParsePara:InvalidInput',...
            'Matrix %s_%d has to be numeric', C_name, i);
    end
end

% check real
for i=1:K
    if (~isreal(C{i}))
        error('MATLAB:Q_DT_SMK_ParsePara:InvalidInput',...
            '%s_%d has to be a real matrix',C_name,i);
    end    
end

% check dimension
for i=1:K
    if (mod( size(C{i},2),size(C{i},1)) ~= 0)
        error('MATLAB:Q_DT_SMK_ParsePara:InvalidInput',...
            '%s_%d is not a set of square matrices',C_name,i);
    end   
end

for i=1:K-1
    if (size(C{i},1) ~= size(C{i+1},1))
        error('MATLAB:Q_DT_SMK_ParsePara:InvalidInput',...
            'The matrices %s_%d and %s_d do not have the same number of rows',C_name,i,C_name,i+1);
    end   
end

% check nonnegativity
for i=1:K
    if (min(min(C{i})) < -10^(-14))
        error('MATLAB:Q_DT_SMK_ParsePara:InvalidInput',...
            'The matrix %s_%d contains negative data',C_name,i);
    end    
end    

% check stochasticity
if processType == 0
    Csum = cell(1,K);
    m = size(C{1},1);
    tmax = zeros(1,K);
    for i = 1:K
        tmax(i) = size(C{i},2)/m;
        Csum{i} = reshape(sum(reshape(C{i}, m*m, tmax(i)), 2), m, m);
        if (max(sum(Csum{i},2)) > 1+10^(-14))||(min(sum(Csum{i},2)) < 1-10^(-14))
            error('MATLAB:Q_DT_SMK_ParsePara:InvalidInput',...
                'The matrix %s_%d(1)+%s_%d(2)...+%s_%d(%d) has to be stochastic', C_name,i,C_name,i,C_name,i,tmax(i));
        end 
    end
    % check that transition matrices Csum are equal for all types
    for i = 1:K-1
        if (max(max(abs(Csum{i+1}-Csum{i}))) > 10^(-14))
        error('MATLAB:Q_DT_SMK_ParsePara:InvalidInput',...
            'The matrices %s_%d(1)+%s_%d(2)...+%s_%d(%d) and %s_%d(1)+%s_%d(2)...+%s_%d(%d)',...
                        'have to be identical ', C_name,i,C_name,i,C_name,i,tmax(i),C_name,i+1,C_name,i+1,C_name,i+1,tmax(i+1));
        end
    end    
else
    %check stochasticity of the sum of the arrival matrices
    m = size(C{1},1);
    Csum = sum(C{1}, 2);
    for i = 2:K
        Csum = Csum + sum(C{i}, 2);
    end
    if (max(Csum) > 1+10^(-14))||(min(Csum) < 1-10^(-14))
      error('MATLAB:Q_DT_SMK_ParsePara:InvalidInput',...
        'The transition matrix of the embedded Markov chain of the semi-Markov process has to be stochastic');
    end 
end