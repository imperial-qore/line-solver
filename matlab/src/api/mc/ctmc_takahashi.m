function [p,p_1,pcourt,Qperm,eps,epsMAX]=ctmc_takahashi(Q,MS,numSteps)
% CTMC_TAKAHASHI - Takahashi's aggregation-disaggregation method
% [p,p_1,pcourt,Qperm,eps,epsMAX] = CTMC_TAKAHASHI(Q,MS,numSteps)
% -- Input
% Q       : infinitesimal generator matrix
% MS      : cell array where MS{i} is the set of rows of Q in macrostate i
% numSteps: number of iterative steps
% -- Output
% p       : estimated steady-state probability vector
% p_1     :
% pcourt  : steady-state probability vector estimated by ctmc_courtois
% Qperm   : permuted Q matrix w.r.t MS as returned by ctmc_courtois
% eps     : NCD index as returned by ctmc_courtois
% epsMAX  : max acceptable value for eps (otherwise Q is not NCD)
% -- Remarks
% * The initial approximate solutions is obtained by calling CTMC_COURTOIS(Q,MS)
% * No convergence stop criterion is currently implented

%% INIT
nMacroStates = size(MS,1); % Number of macro-states
nStates=size(Q,1);
%% START FROM COURTOIS DECOMPOSITION SOLUTION
[pcourt,Qperm,Qdec,eps,epsMAX,P,B,C,q]=ctmc_courtois(Q,MS);
P=ctmc_randomization(Q);
%% STEP 0
pn=pcourt;
%% MAIN LOOP
for n=1:numSteps
pn_1=pn;
p_1=pn_1;
%% AGGREGATION STEP
G=zeros(nMacroStates,nMacroStates);
procRows=0; %processed rows
for I = 1:nMacroStates  % for each source macro-state
    S=sum(pn_1((procRows+1):(procRows+length(MS{I}))));
    procCols=0; %processed cols
    for J = 1:nMacroStates  % for dest macro-state
        if I~=J
            for i=1:length(MS{I})
                for j=1:length(MS{J})
                    if S>1e-14
                    G(I,J)=G(I,J)+P(procRows+i,procCols+j)*pn_1(procRows+i)/S;
                    end
                end
            end
        end
        procCols = procCols + length(MS{J});
    end
    procRows = procRows + length(MS{I});
end
for i = 1:nMacroStates  % for each source macro-state
    G(i,i)=1-sum(G(i,:));
end
gamma=dtmc_solve(G); %compute macroprobabilities

%% DISAGGREGATION STEP
procRows=0; %processed rows
GI=zeros(nMacroStates,nStates);
for I = 1:nMacroStates  % for each source macro-state
    S=sum(pn_1((procRows+1):(procRows+length(MS{I}))));
    for j=1:nStates
        GI(I,j)=GI(I,j)+sum(P((procRows+1):(procRows+length(MS{I})),j).*pn_1((procRows+1):(procRows+length(MS{I})))')/S;
    end    
    procRows = procRows + length(MS{I});
end
procRows=0; %processed rows
for I = 1:nMacroStates  % for each source macro-state
    A=eye(length(MS{I}),length(MS{I}));
    b=zeros(length(MS{I}),1);
    for i=1:length(MS{I})
        for j=1:length(MS{I})
        A(i,j)=A(i,j)-P(procRows+j,procRows+i);
        end
        for K=1:nMacroStates
            if K~=I
                b(i)=b(i)+gamma(K)*GI(K,procRows+i);
            end
        end
    end        
    pn((procRows+1):(procRows + length(MS{I})))=A\b;  
    procRows = procRows + length(MS{I});
end
%% END LOOP
pn=pn/sum(pn);
end
%% OUTPUT
p=pn(:)';
end