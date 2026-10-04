function [sn, options, isFJ] = fjAugment(self, sn, options)
% [SN, OPTIONS, ISFJ] = FJAUGMENT(SELF, SN, OPTIONS)
%
% Apply the tag augmentation a fork-join model needs before its chain can be
% enumerated, and switch the generator to reachability mode. Returns the struct
% unchanged, and ISFJ false, for a model without a Fork or a Join.
%
% A fork firing does not conserve the per-chain population, so the population
% lattice State.spaceGenerator walks cannot enumerate a fork-join model: it
% produces an EMPTY local space for the Join and hence a 0x0 chain. runAnalyzer
% has always augmented first (see its isFJ block); getGenerator and
% getStateSpace did not, so they returned that empty chain instead of an answer
% or an error. See BUGS.md BUG-88.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

isFJ = any(sn.nodetype == NodeType.Fork) || any(sn.nodetype == NodeType.Join);
if ~isFJ
    return
end
[~, fjsn] = ModelAdapter.fjtag(self.model);
sn = fjsn;
options.config.state_space_gen = 'reachable';
end
