function QBD_RAP_ParsePara(A0,A1,A2)
%QBD_RAP_ParsePara checks the validity of the input matrices A0, A1 and A2
%for a QBD with RAP components. The evaluated conditions are necessary but not
%sufficient.

% check numeric 
if (~isnumeric(A0))
    error('MATLAB:QBD_RAP_ParsePara:InvalidInput',...
        'A0 has to be numeric');
end    
if (~isnumeric(A1))
    error('MATLAB:QBD_RAP_ParsePara:InvalidInput',...
        'A1 has to be numeric');
end    
if (~isnumeric(A2))
    error('MATLAB:QBD_RAP_ParsePara:InvalidInput',...
        'A2 has to be numeric');
end    

% check real
if (~isreal(A0))
    error('MATLAB:QBD_RAP_ParsePara:InvalidInput',...
        'A0 has to be a real matrix');
end    
if (~isreal(A1))
    error('MATLAB:QBD_RAP_ParsePara:InvalidInput',...
        'A1 has to be a real matrix');
end    
if (~isreal(A2))
    error('MATLAB:QBD_RAP_ParsePara:InvalidInput',...
        'A2 has to be a real matrix');
end    

% check dimension
if (size(A0,1) ~= size(A0,2))
    error('MATLAB:QBD_RAP_ParsePara:InvalidInput',...
        'A0 is not a square matrix');
end   
if (size(A1,1) ~= size(A1,2))
    error('MATLAB:QBD_RAP_ParsePara:InvalidInput',...
        'A1 is not a square matrix');
end   
if (size(A2,1) ~= size(A2,2))
    error('MATLAB:QBD_RAP_ParsePara:InvalidInput',...
        'A2 is not a square matrix');
end   
if (size(A0,1) ~= size(A1,1))
    error('MATLAB:QBD_RAP_ParsePara:InvalidInput',...
        'The matrices A0 and A1 do not have the same dimension');
end   
if (size(A0,1) ~= size(A2,1))
    error('MATLAB:QBD_RAP_ParsePara:InvalidInput',...
        'The matrices A0 and A2 do not have the same dimension');
end   

% check zero row sum
if (max(sum(A0+A1+A2,2)) > 10^(-14)) || (min(sum(A0+A1+A2,2)) < -10^(-14))
    error('MATLAB:QBD_RAP_ParsePara:InvalidInput',...
        'The matrix A0+A1+A2 must have zero row sum');
end 

%check dominant eigenvalue 
if (max(real(eig(A1))) > -10^(-14))
    error('MATLAB:QBD_RAP_ParsePara:InvalidInput',...
        'The dominant eigenvalue of the matrix A1 must have negative real part');
end 
if (max(real(eig(A0+A1+A2))) > 10^(-14)) || (max(real(eig(A0+A1+A2))) < -10^(-14))
    error('MATLAB:QBD_RAP_ParsePara:InvalidInput',...
        'The dominant eigenvalue of the matrix A0+A1+A2 must have zero real part');
end 


