function varargout = getAvgSysTable(self,varargin)
% [AVGSYSCHAINTABLE, CT,XT] = GETAVGSYSTABLE(SELF,R,T)
% The result recorder captures the returned table together with the solver
% that produced it, so cross-codebase parity is asserted against the values a
% solver RETURNED rather than the text it printed. Off unless a run asked for
% it (LineResultRecorder.enable), and then it costs one appdata lookup here.
% The wrapper exists so that recording happens on EVERY exit path, including
% the early returns inside the implementation below.
[scope, scopeGuard] = LineResultRecorder.enter(); %#ok<ASGLU>
[varargout{1:max(nargout,1)}] = getAvgSysTable_impl(self,varargin{:});
LineResultRecorder.capture(scope, self, 'sys', varargout{1});
end

function [AvgSysChainTable, CT,XT] = getAvgSysTable_impl(self,R,T)
% GETAVGSYSTABLE_IMPL Implementation of GETAVGSYSTABLE; see the wrapper above.

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
