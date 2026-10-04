function SS = spaceClosedMultiCS(M, N, chains, caps)
% SS = SPACECLOSEDMULTICS(M, N, CHAINS)
% SS = SPACECLOSEDMULTICS(M, N, CHAINS, CAPS)
%
% CAPS, when given, is an MxK per (node, class) capacity handed through to
% spaceClosedMulti so the per-class distribution is pruned at the branch.

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% State space for closed multiclass CQN with class-switching

if nargin < 4
    caps = [];
end
C = size(chains,1);
K = size(chains,2);
chainInitPos = cell(1,C);
inchain = cell(1,C);
for c=1:C
    inchain{c}=find(chains(c,:));
    chainInitPos{c} = multichoose(length(inchain{c}),sum(N(inchain{c})));
end
SS = [];
chainInitPosLen = cellfun(@(c) size(c,1),chainInitPos)-1;
v = pprod(chainInitPosLen);
while v>=0
    % SCATTER each chain's per-class populations into their GLOBAL class
    % positions. Appending them chain by chain instead indexes subN by position
    % within the concatenation of the chains, which agrees with the global class
    % order only when every chain owns a CONTIGUOUS block of classes in order.
    % With one class per item (a cache network) chain 1 owns classes 1-4 and 9,
    % so the appended form permuted the populations onto the wrong classes: it
    % emitted marginals holding two jobs of a population-1 chain and never
    % emitted the one holding a job of each chain, which is the declared initial
    % state.
    subN = zeros(1,K);
    for c=1:C
        subN(inchain{c}) = chainInitPos{c}(v(c)+1,:);
    end
    SS = [SS; State.spaceClosedMulti(M, subN, caps)];
    v = pprod(v,chainInitPosLen);
end
end
