function Q_CT_PH_ParsePara(alpha,alpha_name,A,A_name)
%Q_CT_PH_ParsePara checks the validity of the input parameters alpha, A as a
%representation of a Continuous-Time PH distribution. 

% check numeric 
if (~isnumeric(alpha))
    error('MATLAB:Q_CT_PH_ParsePara:InvalidInput',...
        '%s has to be numeric', alpha_name);
end    
if (~isnumeric(A))
    error('MATLAB:Q_CT_PH_ParsePara:InvalidInput',...
        '%s has to be numeric', A_name);
end    

% check real
if (~isreal(alpha))
    error('MATLAB:Q_CT_PH_ParsePara:InvalidInput',...
        '%s has to be a real vector', alpha_name);
end    
if (~isreal(A))
    error('MATLAB:Q_CT_PH_ParsePara:InvalidInput',...
        '%s has to be a real matrix', A_name);
end    

% check dimension
if (size(alpha,1) ~= 1)
    error('MATLAB:Q_CT_PH_ParsePara:InvalidInput',...
        '%s is not a row vector', alpha_name);
end   
if (size(A,1) ~= size(A,2))
    error('MATLAB:Q_CT_PH_ParsePara:InvalidInput',...
        '%s is not a square matrix', A_name);
end   
if (size(alpha,2) ~= size(A,1))
    error('MATLAB:Q_CT_PH_ParsePara:InvalidInput',...
        'The vector %s and the matrix %s do not have the same number of columns', alpha_name, A_name);
end   

% check nonnegativity of vector alpha
if (min(alpha) < -10^(-14))
    error('MATLAB:Q_CT_PH_ParsePara:InvalidInput',...
        'The vector %s contains negative data', alpha_name);
end 
% check negativity of the diagonal entries of A
if (max(diag(A)) > -10^(-14))
    error('MATLAB:Q_CT_PH_ParsePara:InvalidInput',...
        'Some diagonal entries of matrix %s are not negative', A_name);
end 
% check nonnegativity of the off-diagonal entries of A
if (min(min(A-diag(diag(A)) )) < -10^(-14))
    error('MATLAB:Q_CT_PH_ParsePara:InvalidInput',...
        'Some off-diagonal entries of matrix %s are negative', A_name);
end

% check (sub)stochasticity
if (sum(alpha,2) > 1+10^(-14))
    error('MATLAB:Q_CT_PH_ParsePara:InvalidInput',...
        'The vector %s has to be (sub)stochastic', alpha_name);
end 

% check negative row sum
if (max(sum(A,2)) > 10^(-14))
    error('MATLAB:Q_CT_MAP_ParsePara:InvalidInput',...
        'The matrix %s must have negative row sum', A_name);
end