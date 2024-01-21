function [QN,UN,RN,TN] = getAvgChain(self,~,~,~,~)
% [QN,UN,RN,TN] = GETAVGCHAIN(SELF,~,~,~,~)

% Return average station metrics aggregated by chain
%
% Copyright (c) 2012-2024, Imperial College London
% All rights reserved.
QN = self.getAvgQLenChain;
UN = self.getAvgUtilChain;
RN = self.getAvgRespTChain;
TN = self.getAvgTputChain;
end
