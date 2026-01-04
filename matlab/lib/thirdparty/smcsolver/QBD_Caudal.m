function eta=QBD_CAUDAL(A0,A1,A2,varargin)
%QBD_CAUDAL Computes the Spectral Radius of R 
%
%   eta=QBD_CAUDAL(A0,A1,A2) computes the dominant eigenvalue of the
%   matrix R, the smallest nonnegative solution to R= A2 + R A1 + R^2 A0,
%   in case the QBD is recurrent 
%
%   Optional Parameters:
%   
%       Dual: When set to 1, the dominant eigenvalue of the Ramaswami
%             dual is returned. The input QBD must be transient.

OptionNames=['Dual        '];
OptionTypes=['numeric'];
OptionValues=[];
 
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.Dual=0;

% Parse Optional Parameters
options=ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);

if (options.Dual == 1) 
    % the dominant eigenvalue of the Ramaswami dual is
    % identical to the dominant eigenvalue of the process
    % where A2 and A0 are switched.
    A2old=A2;
    A2=A0;
    A0=A2old;
end    

eta_min=0;
eta_max=1;
eta=1/2;
while (eta_max - eta_min > 10^(-15))
    new_eta=max(eig(A2+A1*eta+A0*eta^2));
    if (new_eta > eta)
        eta_min=eta;
    else 
        eta_max=eta;
    end
    eta=(eta_min+eta_max)/2;
end

