function [AvgSysChainTable, CT,XT] = getAvgSysTable(self,R,T)
% [AVGSYSCHAINTABLE, CT,XT] = GETAVGSYSTABLE(SELF,R,T)

% Return table of average system metrics
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if GlobalConstants.DummyMode
    [AvgSysChainTable, CT, XT] = deal(Table());
    AvgSysChainTable = IndexedTable(AvgSysChainTable);
    return
end

if nargin==1
    R = self.getAvgRespTHandles;
    T = self.getAvgTputHandles;
end

if nargin == 2
    if iscell(R) && ~isempty(R)
        param = R;
        [R, T] = deal(param{1:2});   
        % case where varargin is passed as input
    elseif iscell(R) && isempty(R)
        R = self.getAvgRespTHandles;
        T = self.getAvgTputHandles;
    end
end

[SysRespT, SysTput] = getAvgSys(self, R, T);
SysRespT=SysRespT';
SysTput=SysTput';
ChainObj = self.model.getChains();
Chain = cellfun(@(c) c.name,ChainObj,'UniformOutput',false)';
JobClasses = cell(0,1);
for c=1:length(Chain)    
    JobClasses(c,1) = {label(ChainObj{c}.classnames)};
end
Chain = label(Chain);
CT = Table(Chain, JobClasses, SysRespT);
XT = Table(Chain, JobClasses, SysTput);
AvgSysChainTable = Table(Chain, JobClasses, SysRespT, SysTput);
AvgSysChainTable = IndexedTable(AvgSysChainTable);
end
