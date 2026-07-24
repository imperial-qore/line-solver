function [ql]=Q_RAP_RAP_1(C0,C1,D0,D1,varargin)
%   [ql]=Q_RAP_RAP_1(C0,C1,D0,D1) computes the Queue-length distribution 
%   of a RAP(C0,C1)/RAP(D0,D1)/1/FCFS queue    
%   
%   INPUT PARAMETERS:
%   * RAP(C0,C1) arrival process (with mA states)
%     
%   * RAP(D0,D1) service process (with mS states)
%     
%   RETURN VALUES:
%   * Queue length distribution, 
%     ql(i) = Prob[(i-1) customers in the queue]
%
%   OPTIONAL PARAMETERS:
%       Mode: The underlying function to compute the R matrix of the 
%           underlying QBD can be selected using the following 
%           parameter values (default: 'CR')
%               'CR' : Cyclic Reduction [Bini, Meini]
%               'FI' : Functional Iterations [Neuts]
%               'IS' : Invariant Subspace [Akar, Sohraby]
%               'LR' : Logaritmic Reduction [Latouche, Ramaswami]
%               'NI' : Newton Iteration
%
%       MaxNumComp: Maximum number of components for the vectors containig
%           the performance measure.
%       
%       Verbose: When set to 1, the progress of the computation is printed
%           (default:0).
%
%       Optfname: Optional parameters for the underlying function fname.
%           These parameters are included in a cell with one entry holding
%           the name of the parameter and the next entry the parameter
%           value. In this function, fname can be equal to:
%               'QBD_CR' : Options for Cyclic Reduction [Bini, Meini]
%               'QBD_FI' : Options for Functional Iterations [Neuts]
%               'QBD_IS' : Options for Invariant Subspace [Akar, Sohraby]
%               'QBD_LR' : Options for Logaritmic Reduction [Latouche, Ramaswami]
%               'QBD_NI' : Options for Newton Iteration
%
%   USES: QBD Solver and QBD_pi of the SMCSolver tool


OptionNames=[
             'Mode              ';
             'MaxNumComp        ';
             'Verbose           ';
             'OptQBD_CR         ';
             'OptQBD_FI         ';
             'OptQBD_IS         ';
             'OptQBD_LR         ';
             'OptQBD_NI         '];
         
OptionTypes=[
             'char   ';
             'numeric';
             'numeric';
             'cell   '; 
             'cell   '; 
             'cell   '; 
             'cell   '; 
             'cell   '];

OptionValues{1}=['CR';
                 'FI';
                 'IS';
                 'LR';
                 'NI'];
 
options=[];
for i=1:size(OptionNames,1)
    options.(deblank(OptionNames(i,:)))=[];
end    

% Default settings
options.Mode='CR';
options.MaxNumComp = 50;
options.Verbose = 0;
options.OptQBD_CR=cell(0);
options.OptQBD_FI=cell(0);
options.OptQBD_IS=cell(0);
options.OptQBD_LR=cell(0);
options.OptQBD_NI=cell(0);

% Parse Parameters
Q_RAP_ParsePara(C0,'C0',C1,'C1')
Q_RAP_ParsePara(D0,'D0',D1,'D1')

% Parse Optional Parameters
options=Q_ParseOptPara(options,OptionNames,OptionTypes,OptionValues,varargin);
% Check for unused parameter
Q_CheckUnusedParaQBD(options)

% Determine constants
mA=size(C0,1);
mS=size(D0,1);

% Test the load of the queue
pi_a = stat(C0+C1+eye(mA));
pi_s = stat(D0+D1+eye(mS));
lambda = sum(pi_a*C1);
mu = sum(pi_s*D1);

load=lambda/mu
if load >= 1
    error('MATLAB:Q_RAP_RAP_1:LoadExceedsOne',...
                        'The load %d of the system exceeds one',load);
end    

% Compute QBD blocks A0, A1 and A2
A0=kron(eye(mA),D1);
A1=kron(C0,eye(mS))+kron(eye(mA),D0);
A2=kron(C1,eye(mS));

B0=A0;
B1=kron(C0,eye(mS));

if  (strfind(options.Mode,'FI')>0)
    [G,R]=QBD_FI(A0,A1,A2,options.OptQBD_FI{:},'RAPComp',1);
elseif (strfind(options.Mode,'LR')>0)
    [G,R]=QBD_LR(A0,A1,A2,options.OptQBD_LR{:},'RAPComp',1);
elseif (strfind(options.Mode,'IS')>0)
    [G,R]=QBD_IS(A0,A1,A2,options.OptQBD_IS{:},'RAPComp',1);
elseif (strfind(options.Mode,'NI')>0)
    [G,R]=QBD_NI(A0,A1,A2,options.OptQBD_NI{:},'RAPComp',1);
else
    [G,R]=QBD_CR(A0,A1,A2,options.OptQBD_CR{:},'RAPComp',1);
end

pi  = QBD_pi(B0,B1,R, 'MaxNumComp',options.MaxNumComp,'Verbose',options.Verbose,'RAPComp',1);


% compute queue length 
ql  = zeros(1, size(pi, 2)/(mA*mS));
for i = 1:size(pi,2)/(mA*mS)
    ql(i) =sum(pi((i-1)*mA*mS+1:i*mA*mS));
end