function Q_RAP_ParsePara2(C,C_name,D,D_name)
%Q_RAP_ParsePara checks the validity of the input matrices C, D as a
%representation of RAP. The evaluated conditions are necessary but not
%sufficient.

% check numeric 
if (~isnumeric(C))
    error('MATLAB:Q_RAP_ParsePara:InvalidInput',...
        '%s has to be numeric', C_name);
end    
if (~isnumeric(D))
    error('MATLAB:Q_RAP_ParsePara:InvalidInput',...
        '%s has to be numeric', D_name);
end    

% check real
if (~isreal(C))
    error('MATLAB:Q_RAP_ParsePara:InvalidInput',...
        '%s has to be a real matrix', C_name);
end    
if (~isreal(D))
    error('MATLAB:Q_RAP_ParsePara:InvalidInput',...
        '%s has to be a real matrix', D_name);
end    

% check dimension
if (size(C,1) ~= size(C,2))
    error('MATLAB:Q_RAP_ParsePara:InvalidInput',...
        '%s is not a square matrix', C_name);
end   
if (size(D,1) ~= size(D,2))
    error('MATLAB:Q_RAP_ParsePara:InvalidInput',...
        '%s is not a square matrix', D_name);
end   
if (size(C,1) ~= size(D,1))
    error('MATLAB:Q_RAP_ParsePara:InvalidInput',...
        'The matrices %s and %s do not have the same dimension', C_name, D_name);
end   

% check zero row sum
if (max(sum(C+D,2)) > 10^(-14)) || (min(sum(C+D,2)) < -10^(-14))
    error('MATLAB:Q_RAP_ParsePara:InvalidInput',...
        'The matrix %s+%s must have zero row sum', C_name, D_name);
end 

%check dominant eigenvalue 
if (max(real(eig(C))) > -10^(-14))
    error('MATLAB:Q_RAP_ParsePara:InvalidInput',...
        'The dominant eigenvalue of the matrix %s must have negative real part', C_name);
end 
if (max(real(eig(C+D))) > 10^(-14)) || (max(real(eig(C+D))) < -10^(-14))
    error('MATLAB:Q_RAP_ParsePara:InvalidInput',...
        'The dominant eigenvalue of the matrix %s+%s must have zero real part', C_name, D_name);
end 

%check non-negative moments
m = size(C,1);
invC = inv(C);
P = -invC*D;
pi = stat(P);
m1 = -pi*invC*ones(m,1);
if (m1 < -10^(-14))
    error('MATLAB:Q_RAP_ParsePara:InvalidInput',...
        'The first moment of the stationary inter-event distribution of the RAP must be non-negative');
end 
m2 = 2*pi*invC*invC*ones(m,1);
if (m2 < -10^(-14))
    error('MATLAB:Q_RAP_ParsePara:InvalidInput',...
        'The second moment of the stationary inter-event distribution of the RAP must be non-negative');
end 
var = m2 - m1^2;
if (var < -10^(-14))
    error('MATLAB:Q_RAP_ParsePara:InvalidInput',...
        'The variance of the stationary inter-event distribution of the RAP(%s,%s) must be non-negative', C_name, D_name);
end 
m11 = pi*invC*P*invC*ones(m,1);
if (m11 < -10^(-14))
    error('MATLAB:Q_RAP_ParsePara:InvalidInput',...
        'The first joint moment of two succesive inter-events of the stationary RAP(%s,%s) must be non-negative', C_name, D_name);
end 

%check non-negative stationary inter-event density
stdev = sqrt(var);
limSup = m1 + 10*stdev;
quantile = pi*expm(C*limSup)*ones(m,1);
while quantile > 10^-6
    limSup = limSup + stdev;
    quantile = pi*expm(C*limSup)*ones(m,1);
end
numPoints = 200;
for x=1:numPoints
    if( -pi*expm(C*x/numPoints*limSup)*C*ones(m,1) < -10^(-14) )
        error('MATLAB:Q_RAP_ParsePara:InvalidInput',...
        'The stationaty inter-event density function of the RAP(%s,%s) is negative', C_name, D_name);
    end
end

%check non-negative joint density
stdev = sqrt(var);
limSup = m1 + 5*stdev;
quantile = pi*expm(C*limSup)*D*expm(C*limSup)*ones(m,1);
while quantile > 10^-6
    limSup = limSup + stdev;
    quantile = pi*expm(C*limSup)*D*expm(C*limSup)*ones(m,1);
end
numPoints = 100;
for x=1:numPoints
    for y=1:numPoints
        if( pi*expm(C*x/numPoints*limSup)*D*expm(C*y/numPoints*limSup)*D*ones(m,1) < -10^(-14) )
            error('MATLAB:Q_RAP_ParsePara:InvalidInput',...
            'The stationary joint inter-event density function of the RAP(%s,%s) is negative', C_name, D_name);
        end
    end
end

