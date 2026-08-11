function [infGen, eventFilt, syncInfo, stateSpace, nodeStateSpace] = getSymbolicGenerator(self,invertSymbol)
% [INFGEN, EVENTFILT, SYNCINFO, STATESPACE, NODESTATESPACE] = GETSYMBOLICGENERATOR(INVERTSYMBOL)
%
% Symbolic infinitesimal generator, with each event filtration normalized by
% its minimum positive rate and scaled by a symbolic variable x1, ..., xE.
%
% Two backends produce it. With the Symbolic Math Toolbox, INFGEN and the
% entries of EVENTFILT are sym matrices, as they have always been. Without it,
% the same matrices are returned as cell arrays of expression strings built
% through the line-sage-rest service (see SAGE.m), which is what the JAR and
% native Python return as well. The generator is linear in the symbols, so no
% computer algebra is needed to assemble it either way; solving with it does,
% see getSymbolicSolution.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin<2
    invertSymbol = false;
end
if isdeployed
    infGen = [];
    eventFilt = [];
    syncInfo = [];
    stateSpace = [];
    nodeStateSpace = [];
    return
end

useSym = SAGE.hasSymbolicToolbox();
if ~useSym && ~SAGE.isAvailable(SolverCTMC.symbolicBackend(self))
    line_error(mfilename, ['This method requires MATLAB''s Symbolic Toolbox or a ' ...
        'symbolic backend. Start one with: docker run -d -p 8080:8080 imperialqore/line-sage-rest:latest']);
end

[~, F] = getGenerator(self);
[stateSpace, nodeStateSpace] = getStateSpace(self);
n = size(F{1},1);
eventFilt = cell(1, length(F));
if useSym
    infGen = sym(zeros(n));
else
    % Coefficient bookkeeping: the symbolic generator is a sum of numeric
    % matrices scaled by the event symbols, so it is assembled numerically and
    % printed at the end.
    infGenTerms = cell(1, length(F));
    symbols = cell(1, length(F));
end
for e = 1:length(F)
    F{e} = full(F{e});
    minF = min(min(F{e}(F{e}>0)));
    if ~isempty(minF)
        F{e} = F{e} / minF;
        if useSym
            if invertSymbol
                eventFilt{e} = F{e} / sym(['x',num2str(e)],'real');
            else
                eventFilt{e} = F{e} * sym(['x',num2str(e)],'real');
            end
            infGen = infGen + eventFilt{e};
        else
            eventFilt{e} = F{e};
            % ctmc_makeinfgen is linear, so the symbolic generator is the sum
            % of the per-event terms scaled by their symbols.
            infGenTerms{e} = ctmc_makeinfgen(F{e});
            symbols{e} = ['x',num2str(e)];
        end
    end
end
if useSym
    infGen = ctmc_makeinfgen(infGen);
else
    infGen = SolverCTMC.symbolicEntries(infGenTerms, symbols, invertSymbol, n);
end
syncInfo = self.getStruct.sync;
end
