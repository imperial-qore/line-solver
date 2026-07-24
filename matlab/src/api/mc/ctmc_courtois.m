function [p,Qperm,Qdec,eps,epsMAX,P,B,C,q]=ctmc_courtois(Q,MS,q)
% CTMC_COURTOIS - Courtois decomposition
% [p,Qperm,Qdec,eps,epsMAX,P,B,C,q] = CTMC_COURTOIS(Q,MS)
% -- Input
% Q      : infinitesimal generator matrix
% MS     : cell array where MS{i} is the set of rows of Q in macrostate i 
% q      : (optional) randomization coefficient
% -- Output
% p      : approximate steady-state probability vector
% Qperm  : Q reordered according to macrostates
% Qdec   : infinitesimal generator for the macrostates
% P      : probability matrix obtained from Qperm with randomization
% B      : part of P not modelled by decomposition
% eps    : nearly-complete decomposability (NCD) index 
% epsMAX : max acceptable value for eps (otherwise Q is not NCD)
% q      : randomization coefficient

%% INIT
% The randomization coefficient q is derived below, once Qperm is available,
% because it depends on Qperm. It must be a valid uniformization rate
% (q >= max|q_ii|) or P=I+Qperm/q is not stochastic; 1.05*max|Qperm| exceeds
% max|q_ii| strictly, which also keeps P aperiodic. An explicit q (nargin==3)
% overrides the derived value.
v=[];
nMacroStates = size(MS,1); % Number of macro-states
%% REARRANGE INFINITESIMAL GENERATOR ACCORDING TO THE MACROSTATES

for n=1:nMacroStates
    v=[v,MS{n}(:)'];   % force a row: a column MS{n} made v a matrix, and
                       % length(v) then truncated the output p silently
end
Qperm=Q(v,v); % reorder according to the new macro-states
Qdec=Qperm;
procRows=0; %processed rows
for i = 1:size(MS,1)  % for each macro-state
    if procRows >0
        Qdec((procRows + 1):(procRows + length(MS{i})),1:procRows)=0;
    end
    Qdec((procRows + 1):(procRows + length(MS{i})),(procRows + length(MS{i})+1):end)=0;
    procRows = procRows + length(MS{i});
end
% now make each substochastic diagonal block a stochastic matrix
Qdec=ctmc_makeinfgen(Qdec);

%% COMPUTE NCD ERROR INDEX
epsC=Qperm-Qdec;
epsC=0;
C=epsC;
eps=1;

% apply randomization
if nargin==2
    q=(1.05*max(max(abs(Qperm))));
end

P=ctmc_randomization(Qperm,q);
A=P;
procRows=0; %processed rows
for i = 1:nMacroStates % for each macro-state
    if procRows >0
        A((procRows + 1):(procRows + length(MS{i})),1:procRows)=0;
    end
    A((procRows + 1):(procRows + length(MS{i})),(procRows + length(MS{i})+1):end)=0;
    procRows = procRows + length(MS{i});
end
B=P-A;
% see _kb/03-api-layer.md (mc/ additions) for the row-sum-vs-column-sum rationale
eps=max(sum(B,2));
%% COMPUTE epsMAX
% see _kb/03-api-layer.md (mc/ additions) -- eps/epsMAX always computed, not nargout-gated
if true
procRows=0; %processed rows
for i = 1:size(MS,1)  % for each macro-state
    for j=1:length(MS{i})
        pos=j;
        A(procRows+j,procRows+pos)=1-(sum(A(procRows+j,setdiff((procRows+1):(procRows+length(MS{i})),(procRows+pos)))));
    end
    procRows = procRows + length(MS{i});
end
eigMS = zeros(1,size(MS,1));
procRows=0; %processed rows
for i=1:size(MS,1) % for each macro-state
    e=sort(abs(eig(A((procRows + 1):(procRows + length(MS{i})),(procRows + 1):(procRows + length(MS{i}))))));
    if length(e)>1
        eigMS(i)=e(end-1); % take the second largest eigvalues of the block
    else
        eigMS(i)=0; % skip if there is no second eigenvalue
    end
    procRows = procRows + length(MS{i});
end
epsMAX=(1-max(eigMS))/2;
else
    epsMAX=0;   % was epsMax, a typo that never assigned the declared output
end
%% COMPUTE MICROPROBABILITIES
procRows=0; %processed rows
pmicro=zeros(size(Q,1),1);
for i = 1:nMacroStates  % for each macro-state
    Qmicrostate=Qdec((procRows + 1):(procRows + length(MS{i})),(procRows + 1):(procRows + length(MS{i})));
    pmicro((procRows + 1):(procRows + length(MS{i})),1)=ctmc_solve_reducible(Qmicrostate);
    procRows = procRows + length(MS{i});
end

%% COMPUTE MACROPROBABILITIES
G=zeros(nMacroStates,nMacroStates);
procRows=0; %processed rows
for i = 1:nMacroStates  % for each source macro-state
    procCols=0; %processed cols
    for j = 1:nMacroStates  % for dest macro-state
        if i~=j
            for iState=1:length(MS{i})
                G(i,j)=G(i,j)+pmicro(procRows+iState)*sum(P(procRows+iState,(procCols+1):(procCols+length(MS{j}))));
            end            
        end
        procCols = procCols + length(MS{j});
    end
    procRows = procRows + length(MS{i});
end
for i = 1:nMacroStates  % for each source macro-state
    G(i,i)=1-sum(G(i,:));
end
pMacro=dtmc_solve(G);
procRows=0; %processed rows
for i = 1:nMacroStates  % for each source macro-state
    p((procRows+1):(procRows+length(MS{i})))=pMacro(i)*pmicro((procRows+1):(procRows+length(MS{i})));
    procRows = procRows + length(MS{i});
end
%% OUTPUT
for i=1:length(v)
pout(v(i))=p(i);
end
p=pout;

