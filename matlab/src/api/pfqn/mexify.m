%{
%{
 % @brief MATLAB Coder script to generate MEX function for pfqn_bs.
 %
 % This script generates a MEX (MATLAB Executable) version of pfqn_bs
 % for improved performance. It configures the code generation settings
 % and specifies the expected input types and dimensions.
 %
 % See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.
%}
%}
% UNTITLED   Generate static library pfqn_bs
%
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = false;

%% 
ARGS = cell(1,1);
ARGS{1} = cell(7,1);
ARGS{1}{1} = coder.typeof(0,[Inf Inf],[1 1]); %L
ARGS{1}{2} = coder.typeof(0,[1 Inf],[0 1]); %N
ARGS{1}{3} = coder.typeof(0,[1 Inf],[0 1]); %Z
ARGS{1}{4} = coder.typeof(0); %tol
ARGS{1}{5} = coder.typeof(0); %maxiter
ARGS{1}{6} = coder.typeof(0,[Inf Inf],[1 1]); %QN
ARGS{1}{7} = coder.typeof(0,[Inf Inf],[1 1]); %weight
codegen -config cfg pfqn_bs -args ARGS{1}

%% 
%ARGS = cell(1,1);
%ARGS{1} = cell(4,1);
%ARGS{1}{1} = coder.typeof(0,[Inf Inf],[1 1]); %L
%ARGS{1}{2} = coder.typeof(0,[1 Inf],[0 1]); %N
%ARGS{1}{3} = coder.typeof(0,[1 Inf],[0 1]); %Z
%ARGS{1}{4} = coder.typeof(0,[1 Inf],[0 1]); %mi
%codegen -config cfg pfqn_mva -args ARGS{1}
