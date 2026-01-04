% UNTITLED   Generate static library pfqn_bs.
% 
% See also CODER, CODER.CONFIG, CODER.TYPEOF, CODEGEN.

%% Create configuration object of class 'coder.CodeConfig'.
cfg = coder.config('mex','ecoder',false);
cfg.GenerateReport = false;
cfg.ReportPotentialDifferences = false;
cfg.GenCodeOnly = false;
%% 'matchrow'.
ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARGS{1}{1} = coder.typeof(0,[Inf Inf],[1 1]); %A
ARGS{1}{2} = coder.typeof(0,[1 Inf],[0 1]); %B
codegen -config cfg matchrow -args ARGS{1}
%% 'allbut'.
ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARGS{1}{1} = coder.typeof(0,[1 Inf],[0 1]); %y
ARGS{1}{2} = coder.typeof(0); %yps
codegen -config cfg allbut -args ARGS{1}
%% 'at'.
ARGS = cell(1,1);
ARGS{1} = cell(3,1);
ARGS{1}{1} = coder.typeof(0,[Inf Inf],[1 1]); %A
ARGS{1}{2} = coder.typeof(0); %r
ARGS{1}{3} = coder.typeof(0); %i
codegen -config cfg at -args ARGS{1}
%% 'cellat'.
ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARG = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{1} = coder.typeof({ARG}, [1 Inf],[0 1]);
ARGS{1}{2} = coder.typeof(0);
codegen -config cfg cellat -args ARGS{1}
%% 'cellmerge'.
ARGS = cell(1,1);
ARGS{1} = cell(1,1);
ARG = coder.typeof(0,[Inf Inf],[1 1]);
ARGS{1}{1} = coder.typeof({ARG}, [1 Inf],[0 1]);
codegen -config cfg cellmerge -args ARGS{1}
%% 'cellzeros'.
ARGS = cell(1,1);
ARGS{1} = cell(4,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0);
ARGS{1}{4} = coder.typeof(0);
codegen -config cfg cellzeros -args ARGS{1}
%% 'circul'.
% ARGS = cell(1,1);
% ARGS{1} = cell(1,1);
% ARG = coder.typeof(0,[Inf Inf],[1 1]);
% ARGS{1}{1} = coder.typeof(0);
% codegen -config cfg circul -args ARGS{1}
% cannot do, variable size cell array
%% 'findrows'.
ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARGS{1}{1} = coder.typeof(0,[Inf Inf],[1 1]); %M
ARGS{1}{2} = coder.typeof(0,[1 Inf],[0 1]); %r
codegen -config cfg findrows -args ARGS{1}
%% 'factln'.
ARGS = cell(1,1);
ARGS{1} = cell(1,1);
ARGS{1}{1} = coder.typeof(0);
codegen -config cfg factln -args ARGS{1}
%% findstring
% cannot do, variable size cell array
%% 'hashpop'.
ARGS = cell(1,1);
ARGS{1} = cell(4,1);
ARGS{1}{1} = coder.typeof(0,[1 Inf],[0 1]); %r
ARGS{1}{2} = coder.typeof(0,[1 Inf],[0 1]); %r
ARGS{1}{3} = coder.typeof(0); %r
ARGS{1}{4} = coder.typeof(0,[1 Inf],[0 1]); %r
codegen -config cfg hashpop -args ARGS{1}
%% 'krons'.
ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARGS{1}{1} = coder.typeof(0,[Inf Inf],[1 1]); %A
ARGS{1}{2} = coder.typeof(0,[Inf Inf],[1 1]); %B
codegen -config cfg krons -args ARGS{1}
%% 'logmeanexp'.
ARGS = cell(1,1);
ARGS{1} = cell(1,1);
ARGS{1}{1} = coder.typeof(0);
codegen -config cfg logmeanexp -args ARGS{1}
%% 'mape'.
ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0);
codegen -config cfg mape -args ARGS{1}
%% 'maxpe'.
ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0);
codegen -config cfg maxpe -args ARGS{1}
%% 'maxpos'.
ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARGS{1}{1} = coder.typeof(0,[1 Inf],[0 1]); %r
ARGS{1}{2} = coder.typeof(0);
codegen -config cfg maxpos -args ARGS{1}
%% 'multichoosecon'.
% unsupported for find(v) pattern
% ARGS = cell(1,1);
% ARGS{1} = cell(2,1);
% ARGS{1}{1} = coder.typeof(0,[1 Inf],[0 1]); %r
% ARGS{1}{2} = coder.typeof(0);
% codegen -config cfg multichoosecon -args ARGS{1}
%% 'multichoose'.
ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0);
codegen -config cfg multichoose -args ARGS{1}
%% 'multinomialln'.
ARGS = cell(1,1);
ARGS{1} = cell(1,1);
ARGS{1}{1} = coder.typeof(0);
codegen -config cfg multinomialln -args ARGS{1}
%% 'nchoosekln'.
ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0);
codegen -config cfg nchoosekln -args ARGS{1}
%% 'oner'.
ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARGS{1}{1} = coder.typeof(0,[1 Inf],[0 1]); 
ARGS{1}{2} = coder.typeof(0);
codegen -config cfg oner -args ARGS{1}
%% 'pprodcon'.
ARGS = cell(1,1);
ARGS{1} = cell(3,1);
ARGS{1}{1} = coder.typeof(0,[1 Inf],[0 1]); 
ARGS{1}{2} = coder.typeof(0,[1 Inf],[0 1]); 
ARGS{1}{3} = coder.typeof(0,[1 Inf],[0 1]); 
codegen -config cfg pprodcon -args ARGS{1}
%% 'pprod'.
ARGS = cell(1,1);
ARGS{1} = cell(2,1);
ARGS{1}{1} = coder.typeof(0,[1 Inf],[0 1]); 
ARGS{1}{2} = coder.typeof(0,[1 Inf],[0 1]); 
codegen -config cfg pprod -args ARGS{1}
%% 'probchoose'.
ARGS = cell(1,1);
ARGS{1} = cell(1,1);
ARGS{1}{1} = coder.typeof(0,[1 Inf],[0 1]); 
codegen -config cfg probchoose -args ARGS{1}
%% 'softmin'.
ARGS = cell(1,1);
ARGS{1} = cell(3,1);
ARGS{1}{1} = coder.typeof(0);
ARGS{1}{2} = coder.typeof(0);
ARGS{1}{3} = coder.typeof(0);
codegen -config cfg softmin -args ARGS{1}
%% 'sumfinite'.
% ARGS = cell(1,1);
% ARGS{1} = cell(2,1);
% ARG = coder.typeof(0,[Inf Inf],[1 1]);
% ARGS{1}{1} = coder.typeof(0,[1 Inf],[0 1]); %y
% ARGS{1}{2} = coder.typeof(0);
% codegen -config cfg sumfinite -args ARGS{1}
%% 'tget'
% cell arrays
%% 'weaklyconncomp'
% dmperm
