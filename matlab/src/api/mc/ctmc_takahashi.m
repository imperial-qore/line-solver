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
% see _kb/03-api-layer.md (CTMC aggregation-disaggregation) for the randomization rationale
P=ctmc_randomization(Q,1.05*max(max(abs(Q))));
%% STEP 0
pn=pcourt;
% see _kb/03-api-layer.md (CTMC aggregation-disaggregation) for the MS{I} indexing rationale
%% MAIN LOOP
for n=1:numSteps
pn_1=pn;
p_1=pn_1;
%% AGGREGATION STEP
G=zeros(nMacroStates,nMacroStates);
for I = 1:nMacroStates  % for each source macro-state
    idxI=MS{I}(:)';
    S=sum(pn_1(idxI));
    for J = 1:nMacroStates  % for dest macro-state
        if I~=J
            idxJ=MS{J}(:)';
            for i=1:length(idxI)
                for j=1:length(idxJ)
                    if S>1e-14
                    G(I,J)=G(I,J)+P(idxI(i),idxJ(j))*pn_1(idxI(i))/S;
                    end
                end
            end
        end
    end
end
for i = 1:nMacroStates  % for each source macro-state
    G(i,i)=1-sum(G(i,:));
end
gamma=dtmc_solve(G); %compute macroprobabilities

%% DISAGGREGATION STEP
GI=zeros(nMacroStates,nStates);
for I = 1:nMacroStates  % for each source macro-state
    idxI=MS{I}(:)';
    S=sum(pn_1(idxI));
    % The aggregation step above already guards on S>1e-14; dividing here by
    % the same S unguarded put NaN into the iterate whenever a macro-state
    % carried no probability. Zero is the S->0 limit of the contribution.
    if S>1e-14
        for j=1:nStates
            GI(I,j)=GI(I,j)+sum(P(idxI,j).*pn_1(idxI)')/S;
        end
    end
end
for I = 1:nMacroStates  % for each source macro-state
    idxI=MS{I}(:)';
    A=eye(length(idxI),length(idxI));
    b=zeros(length(idxI),1);
    for i=1:length(idxI)
        for j=1:length(idxI)
        A(i,j)=A(i,j)-P(idxI(j),idxI(i));
        end
        for K=1:nMacroStates
            if K~=I
                b(i)=b(i)+gamma(K)*GI(K,idxI(i));
            end
        end
    end
    xI=[];
    if size(A,1) > 6000
        [xg,gflag]=ctmc_gmres(sparse(A),b);
        if gflag==0
            xI=xg;
        end
    end
    if isempty(xI)
        xI=A\b;
    end
    pn(idxI)=xI;
end
%% END LOOP
pn=pn/sum(pn);
end
%% OUTPUT
p=pn(:)';
end