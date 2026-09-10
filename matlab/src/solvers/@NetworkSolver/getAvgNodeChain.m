function [QN,UN,RN,WN,AN,TN] = getAvgNodeChain(self,~,~,~,~)
% [QN,UN,RN,WN,AN,TN] = GETAVGNODECHAIN(SELF,~,~,~,~)

% Return average node metrics aggregated by chain
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
QN = self.getAvgNodeQLenChain;
UN = self.getAvgNodeUtilChain;
RN = self.getAvgNodeRespTChain;
TN = self.getAvgNodeTputChain;
AN = self.getAvgNodeArvRChain;
WN = self.getAvgNodeResidTChain;
end
