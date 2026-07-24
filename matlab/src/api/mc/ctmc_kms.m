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
% see _kb/03-api-layer.md (CTMC aggregation-disaggregation) for the ordering rationale
v=[];
for n=1:nMacroStates
    v=[v,MS{n}];
end
pn=pcourt(v);
p_1=pn(:); % defined even when numSteps==0
%% MAIN LOOP
for n=1:numSteps
    pn_1=pn(:); % canonical column orientation regardless of iteration
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
    zn=zeros(1,nStates); % row: zn*L below requires row orientation
    L=zeros(nStates,nStates);
    D=zeros(nStates,nStates);
    U=zeros(nStates,nStates);
    procRows=0; %processed rows
    for I = 1:nMacroStates  % for each macro-state
        zn((procRows + 1):(procRows + length(MS{I})))=w(I)*z((procRows + 1):(procRows + length(MS{I})))';
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
    % Block Gauss-Seidel sweep: solve pn*(D-U) = zn*L directly. The former
    % 2-block inverse [invA,-invA*B*invC;0,invC] split M at an arbitrary
    % midpoint rather than a macro-block boundary, which is invalid.
    rhs=(zn*L)';
    pn=[];
    if size(M,1) > 6000
        [xg,gflag]=ctmc_gmres(M',rhs);
        if gflag==0
            pn=xg';
        end
    end
    if isempty(pn)
        pn=(M'\rhs)';
    end
    pn=pn/sum(pn);
end
%% OUTPUT
% Map back from permuted (macrostate-major) to original state ordering.
p=zeros(1,nStates);
p(v)=pn;
pback=zeros(nStates,1);
pback(v)=p_1;
p_1=pback;
end