function [QN,UN,RN,WN,AN,TN] = getAvgChain(self,~,~,~,~)
% [QN,UN,RN,WN,TN] = GETAVGCHAIN(SELF,~,~,~,~)

% Return average station metrics aggregated by chain
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
QN = self.getAvgQLenChain;
UN = self.getAvgUtilChain;
RN = self.getAvgRespTChain;
WN = self.getAvgResidTChain;
AN = self.getAvgArvRChain;
TN = self.getAvgTputChain;
end
