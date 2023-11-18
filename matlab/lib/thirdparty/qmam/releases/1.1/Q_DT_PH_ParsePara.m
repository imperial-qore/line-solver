function Q_DT_PH_ParsePara2(alpha,alpha_name,A,A_name)
%Q_DT_PH_ParsePara checks the validity of the input parameters alpha, A as a
%representation of a Discrete-Time PH distribution. 

% check numeric 
if (~isnumeric(alpha))
    error('MATLAB:Q_DT_PH_ParsePara:InvalidInput',...
        '%s has to be numeric', alpha_name);
end    
if (~isnumeric(A))
    error('MATLAB:Q_DT_PH_ParsePara:InvalidInput',...
        '%s has to be numeric', A_name);
end    

% check real
if (~isreal(alpha))
    error('MATLAB:Q_DT_PH_ParsePara:InvalidInput',...
        '%s has to be a real vector', alpha_name);
end    
if (~isreal(A))
    error('MATLAB:Q_DT_PH_ParsePara:InvalidInput',...
        '%s has to be a real matrix', A_name);
end    

% check dimension
if (size(alpha,1) ~= 1)
    error('MATLAB:Q_DT_PH_ParsePara:InvalidInput',...
        '%s is not a row vector', alpha_name);
end   
if (size(A,1) ~= size(A,2))
    error('MATLAB:Q_DT_PH_ParsePara:InvalidInput',...
        '%s is not a square matrix', A_name);
end   
if (size(alpha,2) ~= size(A,1))
    error('MATLAB:Q_DT_PH_ParsePara:InvalidInput',...
        'The vector %s and the matrix %s do not have the same number of columns', alpha_name, A_name);
end   

% check nonnegativity
if (min(alpha) < -10^(-14))
    error('MATLAB:Q_DT_PH_ParsePara:InvalidInput',...
        'The vector %s contains negative data', alpha_name);
end    
if (min(min(A)) < -10^(-14))
    error('MATLAB:Q_DT_PH_ParsePara:InvalidInput',...
        'The matrix %s contains negative data', A_name);
end    

% check (sub)stochasticity
if (sum(alpha,2) > 1+10^(-14))
    error('MATLAB:Q_DT_PH_ParsePara:InvalidInput',...
        'The vector %s has to be (sub)stochastic', alpha_name);
end 
if (max(sum(A,2)) > 1+10^(-14))||(min(sum(A,2)) > 1-10^(-14))
    error('MATLAB:Q_DT_PH_ParsePara:InvalidInput',...
        'The matrix %s has to be substochastic', A_name);
end 
