function [infGen, eventFilt, ev] = getGenerator(self, options)
% [INFGEN, EVENTFILT] = GETGENERATOR()

% [infGen, eventFilt] = getGenerator(self)
% returns the infinitesimal generator of the CTMC and the
% associated filtration for each event

if nargin>1 && islogical(options)
    line_warning(mfilename,'getGenerator(boolean) is now deprecated - remove the boolean argument.\n');
    options = self.getOptions;
elseif nargin<2
    options = self.getOptions;
end

if self.isChainSolver()
    % Chain mode: the generator is the user-supplied one (P-I for a DTMC).
    % There is no event filtration, since transitions carry no event labels.
    if isempty(self.result) || ~isfield(self.result,'infGen')
        self.runAnalyzer(options);
    end
    infGen = self.result.infGen;
    eventFilt = {};
    ev = {};
    return
end

% lang='cpp' takes the generator from line-cli (-a gen plus -a states), never
% from the MATLAB state-space generator: returning MATLAB's chain here would
% label a MATLAB enumeration as a C++ one, and every getter below reads it.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    if isempty(self.result) || ~isfield(self.result,'infGen')
        g = CPPLINE.generator(self.name, self.model, self.options);
        self.result.infGen = g.infGen;
        self.result.space = g.space;
        self.result.spaceAggr = g.spaceAggr;
        % The C++ sends a state as its per-node block widths and not the split
        % matrices; slicing the space here would be this bridge deciding what a
        % node's block is, so the local spaces are left unset rather than guessed.
        self.result.nodeSpace = {};
        self.result.eventFilt = g.eventFilt;
        % Present only when the C++ CLI sent the derived filtrations; left
        % empty otherwise, so the accessors say they are unavailable rather
        % than reporting an all-zero rate as if it had been measured.
        if isfield(g,'auxFilt') && ~isempty(g.auxFilt)
            self.result.auxFilt = g.auxFilt;
        end
        self.result.pi = g.pi;
    end
    infGen = self.result.infGen;
    eventFilt = self.result.eventFilt;
    % The filtration is returned in the CLI's synchronization order. Handing
    % back sn.sync alongside it is only correct if the two enumerate the same
    % events, so the counts are CHECKED: a mismatched pairing would label one
    % engine's filtration matrix with another engine's event descriptor.
    ev = self.getStruct.sync;
    if numel(ev) ~= numel(eventFilt)
        line_error(mfilename, sprintf(['line-cli returned %d event filtration matrices ' ...
            'where the MATLAB struct declares %d synchronizations; the two event ' ...
            'enumerations cannot be paired. Use lang=''matlab'' for getGenerator on ' ...
            'this model.'], numel(eventFilt), numel(ev)));
    end
    return
end

sn = self.getStruct;

% Note: hide_immediate now selectively preserves Cache immediate transitions
% in solver_ctmc.m, so we no longer need to disable it entirely for Cache nodes

% A fork-join model is generated on its tag-augmented copy, as runAnalyzer does.
% Without the augmentation the population lattice yields an empty local space
% for the Join, and the generator below came back 0x0 with no error raised: a
% caller then reads a vacuous chain as if it were the model. See BUGS.md BUG-88.
[sn, options] = self.fjAugment(sn, options);

if isempty(self.result) || ~isfield(self.result,'infGen')
    %line_warning(mfilename,'The model has not been cached. Running SolverCTMC state space generator.');
    [InfGen,StateSpace,StateSpaceAggr,EventFiltration,~,~,sn,AuxFiltration] = solver_ctmc(sn, options);
    if isempty(InfGen)
        line_error(mfilename, ['SolverCTMC generated an EMPTY infinitesimal generator for this model. ' ...
            'A 0x0 generator is a solver limitation, not a valid chain: every caller that loops over its ' ...
            'states would silently do nothing.']);
    end
    self.result.infGen = InfGen;
    self.result.space = StateSpace;
    self.result.spaceAggr = StateSpaceAggr;
    self.result.nodeSpace = sn.space;
    self.result.eventFilt = EventFiltration;
    % The derived filtrations are cached in their OWN field: eventFilt is
    % paired one-to-one with sn.sync (the count is asserted by callers) and
    % is summed as D1 in sample.m, so a START that rides on an existing arc
    % must not join it. See getEventFiltration.
    self.result.auxFilt = AuxFiltration;
end
infGen = self.result.infGen;
eventFilt = self.result.eventFilt;
ev = sn.sync;
end
