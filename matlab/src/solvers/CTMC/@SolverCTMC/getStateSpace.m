function [stateSpace,localStateSpace] = getStateSpace(self, options)
% [STATESPACE, LOCALSTATESPACE] = GETSTATESPACE()
% 
% STATESPACE: MODEL STATE SPACE
% LOCALSTATESPACE: MARGINAL STATE SPACE LOCAL TO EACH NODE

if nargin<2
    options = self.getOptions;
end
sn = self.getStruct;

if isempty(self.result) || ~isfield(self.result,'space')
    [SS,~,qnc] = State.spaceGenerator(sn, options.cutoff, options);
    sn.space = qnc.space;
%     if options.verbose
%         line_printf('\nCTMC state space size: %d states. ',size(SS,1));
%     end
    self.result.space = SS;
    self.result.nodeSpace = qnc.space;
end

stateSpace = self.result.space;

shift = 1;
localStateSpace = cell(1,length(self.result.nodeSpace));
for i=1:length(self.result.nodeSpace)
    localStateSpace{i} = self.result.space(:,shift:(shift+size(self.result.nodeSpace{i},2)-1));
    shift = shift + size(self.result.nodeSpace{i},2);
end
end