function [stateSpace,localStateSpace] = getStateSpace(self, options)
% [STATESPACE, LOCALSTATESPACE] = GETSTATESPACE()
% 
% STATESPACE: MODEL STATE SPACE
% LOCALSTATESPACE: MARGINAL STATE SPACE LOCAL TO EACH NODE

if nargin<2
    options = self.getOptions;
end

if self.isChainSolver()
    % Chain mode: the state space is the one attached to the chain, or the
    % state indices when the user supplied none. There are no local spaces.
    if isempty(self.result) || ~isfield(self.result,'space')
        self.runAnalyzer(options);
    end
    stateSpace = self.result.space;
    localStateSpace = {stateSpace};
    return
end

% lang='cpp' enumerates the state space in line-cli (-a states). The PER-NODE
% COLUMN WIDTHS ride with it (key NodeWidths), so the second output is sliced
% off the flat space exactly as the native tail below slices it: one column
% block per stateful node, every row kept. Returning {stateSpace} instead --
% one cell holding the whole space -- handed a caller indexing it per node the
% GLOBAL space for every node, and any caller reaching past the first cell
% (tut03_repairmen does) died on "Index exceeds the number of array elements".
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    g = CPPLINE.generator(self.name, self.model, self.options);
    self.result.space = g.space;
    self.result.spaceAggr = g.spaceAggr;
    self.result.infGen = g.infGen;
    self.result.eventFilt = g.eventFilt;
    self.result.pi = g.pi;
    stateSpace = g.space;
    if nargout > 1 && isempty(g.nodeWidths)
        line_error(mfilename, ['line-cli returned this state space without NodeWidths, so its ' ...
            'per-node local spaces cannot be sliced from it. Rebuild common/line-cli from this tree.']);
    end
    if ~isempty(g.nodeWidths) && sum(g.nodeWidths) ~= size(stateSpace, 2)
        line_error(mfilename, sprintf(['line-cli returned NodeWidths summing to %d for a state ' ...
            'space %d columns wide; the two describe different enumerations.'], ...
            sum(g.nodeWidths), size(stateSpace, 2)));
    end
    localStateSpace = cell(1, numel(g.nodeWidths));
    shift = 1;
    for i = 1:numel(g.nodeWidths)
        localStateSpace{i} = stateSpace(:, shift:(shift + g.nodeWidths(i) - 1));
        shift = shift + g.nodeWidths(i);
    end
    return
end

sn = self.getStruct;

% A fork-join model is enumerated on its tag-augmented copy, as runAnalyzer
% does: the population lattice below cannot represent a fork firing and would
% return an empty space. The local spaces then belong to the augmented chain,
% which carries one extra stateful node (the Fork). See BUGS.md BUG-88.
[~, ~, isFJ] = self.fjAugment(sn, options);

if isFJ
    % Take the chain from getGenerator rather than enumerating here: the
    % fork-occupied and join-firable states are vanishing and are removed by
    % the stochastic complement in solver_ctmc, so enumerating separately would
    % return a state space with MORE rows than the generator it is meant to
    % label. One producer, one pairing.
    if isempty(self.result) || ~isfield(self.result,'space')
        self.getGenerator(options);
    end
elseif isempty(self.result) || ~isfield(self.result,'space')
    sn = solver_ctmc_pools(sn); % server pools enumerate in the layout solver_ctmc builds the generator on
    [SS,~,qnc] = State.spaceGenerator(sn, options.cutoff, options);
    sn.space = qnc.space;
    if isempty(SS)
        line_error(mfilename, ['SolverCTMC generated an EMPTY state space for this model, so there is no ' ...
            'chain to return. This is a solver limitation, not an empty model: check that every stateful ' ...
            'node admits a local state (a Join outside a fork-join augmentation does not).']);
    end
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