function [p,p_1,pcourt,Qperm,eps,epsMAX]=ctmc_multi(Q,MS,MSS)
% CTMC_MULTI - Multigrid aggregation-disaggregation method
% CTMC_Multigrid is the basic one step implementation of standard method of
% Multigrid. The complete implementation (e.g multi-level disaggregation)
% requires repetitive coarsening using the base method. For illustrative
% purpose, we only use two levels.

% [p,p_1,pcourt,Qperm,eps,epsMAX] = CTMC_multi(Q,MS,MSS)
% -- Input
% Q      : infinitesimal generator matrix
% MS     : cell array where MS{i} is the set of rows of Q in macrostate i
% MSS    : cell array where MSS{i} is the set of rows of G in macromacrostate i

% -- Output
% p      : approximate steady-state probability vector
% p_1    :
% pcourt : steady-state probability vector estimated by ctmc_courtois
% Qperm  : Q reordered according to macrostates
% eps    : nearly-complete decomposability (NCD) index
% epsMAX : max acceptable value for eps (otherwise Q is not NCD)

%% INIT
if nargin==2
    q=1;
end
v=[];
nMacroStates = size(MS,1); % Number of macro-states
%% REARRANGE INFINITESIMAL GENERATOR ACCORDING TO THE MACROSTATES

for n=1:nMacroStates
    v=[v,MS{n}];
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
if nargin==3
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
eps=max(sum(B));
%% COMPUTE epsMAX
% the following subprocedure makes each diagonal block stochastic by
% placing a normalization condition in the diagonal position
if nargout>5
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
    eps=0;
    epsMAX=0;
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
% Note here we work out the solution of macro states using another
% Courtois method, which means we are solving an even coarsener system,
% then work all the way back to find the solution of original problem
pMacro=ctmc_courtois(G,MSS);
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
