function [p,p_1,Qperm,eps,epsMAX,pcourt]=ctmc_kms(Q,MS,numSteps)
% CTMC_KMS - Koury-McAllister-Stewart aggregation-disaggregation method
% [p,p_1,Qperm,eps,epsMAX,pcourt] = CTMC_KMS(Q,MS,numSteps)
% -- Input
% Q       : infinitesimal generator matrix
% MS      : cell array where MS{i} is the set of rows of Q in macrostate i
% numSteps: number of iterative steps
% -- Output
% p       : estimated steady-state probability vector
% p_1
% Qperm   : permuted Q matrix w.r.t MS as returned by ctmc_courtois
% eps     : NCD index as returned by ctmc_courtois
% epsMAX  : max acceptable value for eps (otherwise Q is not NCD)
% pcourt  : steady-state probability vector estimated by ctmc_courtois
% -- Remarks
% * The initial approximate solutions is obtained by calling CTMC_COURTOIS(Q,MS)
% * No convergence stop criterion is currently implented

%% INIT
nMacroStates = size(MS,1); % Number of macro-states
nStates=size(Q,1);
%% START FROM COURTOIS DECOMPOSITION SOLUTION
[pcourt,Qperm,Qdec,eps,epsMAX,P,B,C]=ctmc_courtois(Q,MS);
%% STEP 0
pn=pcourt;
%% MAIN LOOP
for n=1:numSteps
    pn_1=pn';
    p_1=pn_1;
    %% AGGREGATION STEP
    pcondn_1=pn_1;
    procRows=0; %processed rows
    for I = 1:nMacroStates  % for each macro-state
        if sum(pn_1((procRows + 1):(procRows + length(MS{I}))))>0
        pcondn_1((procRows + 1):(procRows + length(MS{I})))=pn_1((procRows + 1):(procRows + length(MS{I})))/sum(pn_1((procRows + 1):(procRows + length(MS{I}))));
        end
        procRows = procRows + length(MS{I});
    end

    G=zeros(nMacroStates,nMacroStates);
    procCols=0; %processed rows
    for I = 1:nMacroStates  % for each source macro-state
        procRows=0; %processed rows
        for J = 1:nMacroStates  % for each dest macro-state
            %I,J
            %size(pcondn_1((procRows + 1):(procRows + length(MS{J})))')
            %size(P((procRows + 1):(procRows + length(MS{J})),(procCols + 1):(procCols + length(MS{I}))))
            %size(ones(length(MS{I}),1))
            G(I,J)=pcondn_1((procRows + 1):(procRows + length(MS{J})))'*P((procRows + 1):(procRows + length(MS{J})),(procCols + 1):(procCols + length(MS{I})))*ones(length(MS{I}),1);
            procRows = procRows + length(MS{J});
        end
        procCols = procCols + length(MS{I});
    end
    w=dtmc_solve(G');

    %% DISAGGREGATION STEP
    z=pcondn_1;
    L=zeros(nStates,nStates);
    D=zeros(nStates,nStates);
    U=zeros(nStates,nStates);
    procRows=0; %processed rows
    for I = 1:nMacroStates  % for each macro-state
        zn((procRows + 1):(procRows + length(MS{I})))=w(I)*z((procRows + 1):(procRows + length(MS{I})));
        procCols=0; %processed rows
        for J = 1:nMacroStates
            if I>J
                L((procRows + 1):(procRows + length(MS{I})),(procCols + 1):(procCols + length(MS{J})))=P((procRows + 1):(procRows + length(MS{I})),(procCols + 1):(procCols + length(MS{J})));
            end
            if I==J
                D((procRows + 1):(procRows + length(MS{I})),(procCols + 1):(procCols + length(MS{J})))=eye(length(MS{I}))-P((procRows + 1):(procRows + length(MS{I})),(procCols + 1):(procCols + length(MS{J})));
            end
            if I<J
                U((procRows + 1):(procRows + length(MS{I})),(procCols + 1):(procCols + length(MS{J})))=P((procRows + 1):(procRows + length(MS{I})),(procCols + 1):(procCols + length(MS{J})));
            end
            procCols = procCols + length(MS{J});
        end
        procRows = procRows + length(MS{I});
    end
    M=(D-U);
    MPL=round((nStates-2)/2);
    A=M(1:(MPL+1),1:(MPL+1)); invA=inv(A);
    B=M(1:(MPL+1),(MPL+2):end);
    C=M((MPL+2):end,(MPL+2):end); invC=inv(C);
    pn=zn*L*[invA,-invA*B*invC;0*A,invC];
end
%% OUTPUT
p=pn;
end