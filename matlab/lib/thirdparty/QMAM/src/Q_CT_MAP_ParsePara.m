function Q_CT_MAP_ParsePara(C,C_name,D,D_name)
%Q_CT_MAP_ParsePara checks the validity of the input matrices C, D as a
%representation of a continuous-time MAP. 

% check numeric 
if (~isnumeric(C))
    error('MATLAB:Q_CT_MAP_ParsePara:InvalidInput',...
        '%d has to be numeric',C_name);
end    
if (~isnumeric(D))
    error('MATLAB:Q_CT_MAP_ParsePara:InvalidInput',...
        '%s has to be numeric',D_name);
end    

% check real
if (~isreal(C))
    error('MATLAB:Q_CT_MAP_ParsePara:InvalidInput',...
        '%s has to be a real matrix',C_name);
end    
if (~isreal(D))
    error('MATLAB:Q_CT_MAP_ParsePara:InvalidInput',...
        '%s has to be a real matrix',D_name);
end    

% check dimension
if (size(C,1) ~= size(C,2))
    error('MATLAB:Q_CT_MAP_ParsePara:InvalidInput',...
        '%s is not a square matrix',C_name);
end   
if (size(D,1) ~= size(D,2))
    error('MATLAB:Q_CT_MAP_ParsePara:InvalidInput',...
        '%s is not a square matrix',D_name);
end   
if (size(C,1) ~= size(D,1))
    error('MATLAB:Q_CT_MAP_ParsePara:InvalidInput',...
        'The matrices %s and %s do not have the same dimension',C_name,D_name);
end   

% check negativity of the diagonal entries of C
if (max(diag(C)) > -10^(-14))
    error('MATLAB:Q_CT_MAP_ParsePara:InvalidInput',...
        'Some diagonal entries of matrix %s are not negative',C_name);
end 
% check nonnegativity of the off-diagonal entries of C
if (min(min(C-diag(diag(C)) )) < -10^(-14))
    error('MATLAB:Q_CT_MAP_ParsePara:InvalidInput',...
        'Some off-diagonal entries of matrix %s are negative',C_name);
end    
% check nonnegativity of D
if (min(min(D)) < -10^(-14))
    error('MATLAB:Q_CT_MAP_ParsePara:InvalidInput',...
        'The matrix %s contains negative data',D_name);
end    
    
% check zero row sum
if (max(sum(C+D,2)) > 10^(-14)) || (min(sum(C+D,2)) < -10^(-14))
    error('MATLAB:Q_CT_MAP_ParsePara:InvalidInput',...
        'The matrix %s+%s must have zero row sum',C_name,D_name);
end 