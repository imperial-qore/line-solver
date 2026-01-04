function Q_CT_MMAPK_ParsePara(D0,D0_name,D,D_name)
%Q_CT_MMAPK_ParsePara checks the validity of the input matrices D0, D as a
%representation of a Continuous-Time MMAP[K]. D{i} holds the mxm matrix D_i.

% check numeric 
if (~isnumeric(D0))
    error('MATLAB:Q_CT_MMAPK_ParsePara:InvalidInput',...
        '%s has to be numeric',D0_name);
end    
K=size(D,2);
for i=1:K
    if (~isnumeric(D{i}))
        error('MATLAB:Q_CT_MMAPK_ParsePara:InvalidInput',...
            '%s_%d has to be numeric',D_name, i);
    end    
end

% check real
if (~isreal(D0))
    error('MATLAB:Q_CT_MMAPK_ParsePara:InvalidInput',...
        '%s has to be a real matrix',D0_name);
end    
for i=1:K
    if (~isreal(D{i}))
        error('MATLAB:Q_CT_MMAPK_ParsePara:InvalidInput',...
            '%s_%d has to be a real matrix',D_name,i);
    end    
end

% check dimension
if (size(D0,1) ~= size(D0,2))
    error('MATLAB:Q_CT_MMAPK_ParsePara:InvalidInput',...
        '%s is not a square matrix',D0_name);
end   
for i=1:K
    if (size(D{i},1) ~= size(D{i},2))
        error('MATLAB:Q_CT_MMAPK_ParsePara:InvalidInput',...
            '%s_%d is not a square matrices',D_name,i);
    end   
end

if (size(D0,1) ~= size(D{1},1))
    error('MATLAB:Q_CT_MMAPK_ParsePara:InvalidInput',...
        'The matrices %s and %s_1 do not have the same dimension',D0_name,D_name);
end   
for i=1:K-1
    if (size(D{i},1) ~= size(D{i+1},1))
        error('MATLAB:Q_CT_MMAPK_ParsePara:InvalidInput',...
            'The matrices %s_%d and %s_d do not have the same dimension',D_name,i,D_name,i+1);
    end   
end

% check negativity of the diagonal entries of D0
if (max(diag(D0)) > 10^(-14))
    error('MATLAB:Q_CT_MMAPK_ParsePara:InvalidInput',...
        'Some diagonal entries of the matrix %s are not negative',D0_name);
end 
% check nonnegativity of the off-diagonal entries of D0
if (min(min(D0-diag(diag(D0)) )) < -10^(-14))
    error('MATLAB:Q_CT_MMAPK_ParsePara:InvalidInput',...
        'Some off-diagonal entries of the matrix %s are negative',D0_name);
end    
% check nonnegativity of matrices D{i}
for i=1:K
    if (min(min(D{i})) < -10^(-14))
        error('MATLAB:Q_CT_MMAPK_ParsePara:InvalidInput',...
            'The matrix %s_%d contains negative data',D_name,i);
    end    
end
% check zero row sum
Dsum = D0;
for i = 1:K
    Dsum = Dsum + D{i};
end
if (max(sum(Dsum,2)) > 10^(-14)) || (min(sum(Dsum,2)) < -10^(-14))
    error('MATLAB:Q_CT_MMAPK_ParsePara:InvalidInput',...
        'The matrix %s+%s_1+...+%s_K must have zero row sum',D0_name,D_name,D_name);
end 