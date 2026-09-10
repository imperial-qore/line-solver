function stateSpaceAggr = getStateSpaceAggr(self)
% STATESPACEAGGR = GETSTATESPACEAGGR()

options = self.getOptions;

if self.isChainSolver()
    % Chain mode: no phases, so the aggregate space is the state space.
    stateSpaceAggr = self.getStateSpace();
    return
end

% lang='cpp' takes the aggregate space off the same -a states payload the
% space and the stationary law come from, so the three are one enumeration.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    g = CPPLINE.generator(self.name, self.model, self.options);
    self.result.space = g.space;
    self.result.spaceAggr = g.spaceAggr;
    self.result.infGen = g.infGen;
    self.result.eventFilt = g.eventFilt;
    self.result.pi = g.pi;
    stateSpaceAggr = g.spaceAggr;
    return
end

% The aggregate space is derivable from the model alone, exactly like the
% space in getStateSpace and the generator in getGenerator, so build it on
% demand rather than returning [] when nothing has been cached.
if options.force || isempty(self.result) || ~isfield(self.result,'spaceAggr')
    sn = self.getStruct;
    [InfGen,StateSpace,StateSpaceAggr,EventFiltration,~,~,sn] = solver_ctmc(sn, options);
    self.result.infGen = InfGen;
    self.result.space = StateSpace;
    self.result.spaceAggr = StateSpaceAggr;
    self.result.nodeSpace = sn.space;
    self.result.eventFilt = EventFiltration;
end
stateSpaceAggr = self.result.spaceAggr;
end