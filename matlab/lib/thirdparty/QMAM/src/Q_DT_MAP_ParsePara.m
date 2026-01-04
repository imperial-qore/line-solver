function Q_DT_MAP_ParsePara(C,C_name,D,D_name)
%Q_DT_MAP_ParsePara checks the validity of the input matrices C, D as a
%representation of a discrete-time MAP. 

% check numeric 
if (~isnumeric(C))
    error('MATLAB:Q_DT_MAP_ParsePara:InvalidInput',...
        '%d has to be numeric',C_name);
end    
if (~isnumeric(D))
    error('MATLAB:Q_DT_MAP_ParsePara:InvalidInput',...
        '%s has to be numeric',D_name);
end    

% check real
if (~isreal(C))
    error('MATLAB:Q_DT_MAP_ParsePara:InvalidInput',...
        '%s has to be a real matrix',C_name);
end    
if (~isreal(D))
    error('MATLAB:Q_DT_MAP_ParsePara:InvalidInput',...
        '%s has to be a real matrix',D_name);
end    

% check dimension
if (size(C,1) ~= size(C,2))
    error('MATLAB:Q_DT_MAP_ParsePara:InvalidInput',...
        '%s is not a square matrix',C_name);
end   
if (size(D,1) ~= size(D,2))
    error('MATLAB:Q_DT_MAP_ParsePara:InvalidInput',...
        '%s is not a square matrix',D_name);
end   
if (size(C,1) ~= size(D,1))
    error('MATLAB:Q_DT_MAP_ParsePara:InvalidInput',...
        'The matrices %s and %s do not have the same dimension',C_name,D_name);
end   

% check nonnegativity
if (min(min(C)) < -10^(-14))
    error('MATLAB:Q_DT_MAP_ParsePara:InvalidInput',...
        'The matrix %s contains negative data',C_name);
end    
if (min(min(D)) < -10^(-14))
    error('MATLAB:Q_DT_MAP_ParsePara:InvalidInput',...
        'The matrix %s contains negative data',D_name);
end    
    

% check (sub)stochasticity
if (max(sum(C,2)) > 1+10^(-14))
    error('MATLAB:Q_DT_MAP_ParsePara:InvalidInput',...
        'The matrix %s has to be substochastic',C_name);
end    
if (max(sum(C+D,2)) > 1+10^(-14))||(min(sum(C+D,2)) < 1-10^(-14))
    error('MATLAB:Q_DT_MAP_ParsePara:InvalidInput',...
        'The matrix %s+%s has to be stochastic',C_name+D_name);
end    
