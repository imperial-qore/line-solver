function entries = citations(self)
% ENTRIES = CITATIONS()
%
% Bibliographic references for the algorithms this solver used, as a struct
% array with fields .ref (the reference) and .covers (which part of the
% solution process it covers), plus .key, the internal bibliography key of
% doc/latex/biblio.bib, kept for programmatic lookup and never displayed.
% Printed as a list when called without an output argument.
%
% The references follow what the run actually did: the method the analyzer
% resolved to (not merely the one requested), the fork-join transformation if
% the model has forks, and the percentile method of the last getPerctRespT
% call. Running an AMVA model with method 'bs' and then asking for a fork-join
% tail with 'forktail' therefore returns both papers.
%
% Attribution in LINE is pull-based: nothing is printed during a solve.
% See also libraries, which reports the bundled third-party dependencies.
%
% Example:
%   solver = SolverMVA(model,'method','bs');
%   solver.getAvg();
%   solver.getPerctRespT([95 99], [], 'forktail');
%   solver.citations()
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

tokens = {};

% the solver family, both as a method name of its own (for the paradigms whose
% reference IS the paradigm) and as a qualifier on the method names, since the
% same name means different algorithms in different solvers ('mva' is
% Reiser-Lavenberg under MVA and Reiser's convolution under NC)
family = '';
switch class(self)
    case 'SolverCTMC'
        family = 'ctmc'; tokens{end+1} = 'ctmc';
    case 'SolverSSA'
        family = 'ssa';  tokens{end+1} = 'ssa';
    case {'SolverFLD','SolverFluid'}
        family = 'fld';  tokens{end+1} = 'fld';
    case 'SolverLN'
        family = 'ln';   tokens{end+1} = 'ln';
    case 'SolverJMT'
        family = 'jmt';  tokens{end+1} = 'jmt';
    case 'SolverLQNS'
        family = 'lqns'; tokens{end+1} = 'lqns';
    case 'SolverENV'
        family = 'env';  tokens{end+1} = 'env';
    case 'SolverMVA'
        family = 'mva';
    case 'SolverNC'
        family = 'nc';
    case 'SolverMAM'
        family = 'mam';
    case 'SolverBA'
        family = 'ba';
end

    function addMethodToken(m)
        m = lower(strtrim(char(m)));
        if isempty(m) || strcmp(m,'default')
            return
        end
        if ~isempty(family)
            tokens{end+1} = [family '.' m];
        else
            tokens{end+1} = m;
        end
    end

% the requested method, and the one the analyzer actually resolved to: the
% reported name is 'default/<actual>' when dispatch chose for the user
options = self.getOptions;
if isfield(options,'method') && ~isempty(options.method)
    addMethodToken(options.method);
end
if self.hasAvgResults() && isfield(self.result,'Avg') && isfield(self.result.Avg,'method') ...
        && ~isempty(self.result.Avg.method)
    % the reported name is 'default/<actual>' when dispatch chose for the user
    parts = strsplit(char(self.result.Avg.method), '/');
    for pi = 1:numel(parts)
        addMethodToken(parts{pi});
    end
end

% THE BETHE ARM CARRIES ITS OBJECTIVE. The polytope, the quadratic reduction
% and the metric readout of 'qrf.bethe' are the QRF paper's; the functional
% minimised over them is the tree-reweighted free entropy, and the weight
% lambda = 1/M is chosen by the spanning-tree polytope condition of Wainwright,
% Jaakkola and Willsky. Both papers are needed to write the run up.
for ti = 1:numel(tokens)
    tk = char(tokens{ti});
    if strcmp(tk,'qrf.bethe') || strcmp(tk,'ba.qrf.bethe') || ...
            strcmp(tk,'qrf.bas.bethe') || strcmp(tk,'ba.qrf.bas.bethe')
        tokens{end+1} = 'qrf.trw'; %#ok<AGROW>
        break
    end
end

% THE ITERATIVE PB(k) AND BJB(k) CARRY THE HIERARCHY THEY EVALUATE. Casale,
% Muntz and Serazzi name the iteration counts and tabulate them, but both
% brackets are produced by the Eager-Sevcik performance bound hierarchy
% recursion (pfqn_pbh), so a run write-up needs that paper beside theirs.
for ti = 1:numel(tokens)
    tk = char(tokens{ti});
    if any(strcmp(tk, {'pbk','pbk.upper','pbk.lower', ...
                       'bjbk','bjbk.upper','bjbk.lower', ...
                       'ba.pbk','ba.pbk.upper','ba.pbk.lower', ...
                       'ba.bjbk','ba.bjbk.upper','ba.bjbk.lower'}))
        tokens{end+1} = 'pbh'; %#ok<AGROW>
        break
    end
end

% THE LOAD-DEPENDENT DIVDIFF ROUTE CARRIES TWO PAPERS. The outer divided
% difference over the class populations is Casale (SIGMETRICS 2017); only the
% single-class kernel it substitutes is the limited load-dependent closed form
% of Casale, Harrison and Ong. A run write-up needs both.
for ti = 1:numel(tokens)
    tk = char(tokens{ti});
    if strcmp(tk,'divdiff.ld') || strcmp(tk,'nc.divdiff.ld')
        tokens{end+1} = 'divdiff'; %#ok<AGROW>
        break
    end
end

% THE RCAT ROUTE CARRIES ITS QBD, and a user writing the run up needs it: every
% isolated component is a quasi-birth-death process over (queue length, phase),
% and its open tail is Neuts' rate matrix R rather than a scalar ratio.
for ti = 1:numel(tokens)
    tk = char(tokens{ti});
    if any(strcmp(tk, {'inap','inapplus','inapinf', ...
                       'ag.inap','ag.inapplus','ag.inapinf', ...
                       'mam.inap','mam.inapplus','mam.inapinf'}))
        tokens{end+1} = 'rcat.qbd'; %#ok<AGROW>
        break
    end
end

% THE DAE ROUTE CARRIES THREE METHODS BESIDE ITS CLOSURE, and a user writing the
% run up needs all three: the Rosenbrock integrator that takes the singular mass
% matrix, the active set that decides which capacity limits bind, and -- when a
% finite horizon was asked for with a cap present -- the event location that cuts
% the trajectory into segments.
isDae = false;
for ti = 1:numel(tokens)
    tk = char(tokens{ti});
    if strcmp(tk,'dae') || (numel(tk) > 4 && strcmp(tk(end-3:end), '.dae'))
        isDae = true;
        break
    end
end
if isDae
    tokens{end+1} = 'dae.integrator';
    if isa(self.model,'Network')
        snc = self.model.getStruct;
        % a cap the population cannot reach is not a cap: REFRESHCAPACITY derives a
        % FINITE classcap at every station of every closed model, so finiteness
        % alone would report the active set for models that have no constraint
        hasCap = isfield(snc,'nregions') && ~isempty(snc.nregions) && snc.nregions > 0;
        if ~hasCap && isfield(snc,'cap')
            hasCap = any(isfinite(snc.cap(:)) & snc.cap(:) < sum(snc.njobs));
        end
        if ~hasCap && isfield(snc,'classcap')
            for rr = 1:min(size(snc.classcap,2), numel(snc.njobs))
                if any(isfinite(snc.classcap(:,rr)) & snc.classcap(:,rr) < snc.njobs(rr))
                    hasCap = true;
                    break
                end
            end
        end
        if hasCap
            tokens{end+1} = 'dae.activeset';
            if isfield(options,'timespan') && numel(options.timespan) > 1 && isfinite(options.timespan(2))
                tokens{end+1} = 'dae.events';
            end
        end
    end
end

% THE LDQBD METHOD CARRIES A SECOND CONSTRUCTION when the queue is a multiserver
% with phase-type service: the level's inner coordinate is then the MULTISET of
% the phases the busy servers sit in, which is Asmussen and Moller's state space
% rather than Phung-Duc's recursion. Only that shape uses it -- exponential
% service, or a single server, needs no configuration coordinate at all.
isLdqbd = false;
for ti = 1:numel(tokens)
    tk = char(tokens{ti});
    if strcmp(tk,'ldqbd') || (numel(tk) > 6 && strcmp(tk(end-5:end), '.ldqbd'))
        isLdqbd = true;
        break
    end
end
if isLdqbd && isa(self.model,'Network')
    snq = self.model.getStruct;
    % a Delay carries nservers = Inf and the open Source is exponential by the
    % method's own guard, so a finite c > 1 with a non-EXP process is the queue
    isMultiserver = any(isfinite(snq.nservers(:)) & snq.nservers(:) > 1);
    isPHservice = any(any(snq.procid ~= ProcessType.EXP & isfinite(snq.rates) & snq.rates > 0));
    if isMultiserver && isPHservice
        tokens{end+1} = 'ldqbd_mphc';
    end
end

% THE TBI ROUTE CARRIES ITS RELAXATION SCHEME: the cell decomposition is one
% method and the sweep that reconciles the cells is another, so a user writing
% the run up needs Lelarasmee's waveform relaxation beside the TBI paper.
isTbi = false;
for ti = 1:numel(tokens)
    tk = char(tokens{ti});
    if strcmp(tk,'tbi') || (numel(tk) > 4 && strcmp(tk(end-3:end), '.tbi'))
        isTbi = true;
        break
    end
end
if isTbi
    tokens{end+1} = 'tbi.relaxation';
end

% THE COUPLED LAYERED TRANSIENT IS THE SAME RELAXATION over the LQN ensemble.
% It only runs when a finite horizon was asked for: without one getTranAvg
% falls back to the decoupled, frozen-demand transient and no sweep happens.
if strcmp(family,'ln') && isfield(options,'timespan') && numel(options.timespan) >= 2 ...
        && all(isfinite(options.timespan))
    lnTran = 'coupled'; % SolverLN.getTranAvg default
    if isfield(options,'config') && isfield(options.config,'ln_transient') ...
            && ~isempty(options.config.ln_transient)
        lnTran = char(options.config.ln_transient);
    end
    if strcmpi(lnTran,'coupled')
        tokens{end+1} = 'ln.transient.coupled';
    end
end

% the layering strategy of a layered solve, when it is not the default one
if strcmp(family,'ln') && isfield(options,'config') && isfield(options.config,'layering') ...
        && any(strcmpi(options.config.layering, {'flat','squashed'}))
    tokens{end+1} = 'ln.flat';
end

% the encoding the layered solve actually BUILT, which the 'srvn' alias resolves
% at layer build time and which options.method therefore does not name
if strcmp(family,'ln') && isprop(self,'lnmethod') && ~isempty(self.lnmethod)
    addMethodToken(self.lnmethod);
end

% the fork-join transformation, when the model has one
if isa(self.model,'Network') && self.model.hasFork
    fjm = 'mmt';
    if isfield(options,'config') && isfield(options.config,'fork_join') && ~isempty(options.config.fork_join)
        fjm = options.config.fork_join;
        if any(strcmpi(fjm, {'default','fjt'}))
            fjm = 'mmt';
        elseif strcmpi(fjm, 'heidelberger-trivedi')
            fjm = 'ht';
        end
    end
    tokens{end+1} = fjm;
    % a quorum join is a second method on top of the transformation: the
    % synchronisation delay it charges is an order statistic, not a maximum
    if sn_has_quorum_join(self.model.getStruct)
        tokens{end+1} = 'quorum';
    end
end

% the random-environment image of the MAP/MMPP processes, when the run was
% dispatched through it rather than solving the model natively
if self.hasAvgResults() && isfield(self.result,'Avg') && isfield(self.result.Avg,'method') ...
        && contains(char(self.result.Avg.method), 'env.')
    tokens{end+1} = 'env';
    tokens{end+1} = 'map2renv';
    if contains(char(self.result.Avg.method), 'env.dec')
        tokens{end+1} = 'env.dec';
    else
        tokens{end+1} = 'env.avg';
    end
end

% the percentile method of the last getPerctRespT call, if any
if ~isempty(self.lastPerctMethod)
    tokens{end+1} = self.lastPerctMethod;
end

% the permanent engine of the last getProbSysMarg call, if any. The identity
% behind the metric is Ryser's expansion in every case, so 'perm' is reported
% alongside the estimator that evaluated it.
if ~isempty(self.lastPermEngine)
    tokens{end+1} = 'perm';
    if ~strcmpi(self.lastPermEngine,'exact')
        tokens{end+1} = ['perm.' lower(self.lastPermEngine)];
    end
end

entries = line_citations(tokens);

if nargout == 0
    if isempty(entries)
        line_printf('No algorithm references recorded for this run.\n');
    else
        % The bibliography key (.key) is kept on the struct for programmatic
        % lookup in doc/latex/biblio.bib, but it is internal shorthand and is
        % not shown: a user reads the reference, not 'Sch79'.
        for i = 1:numel(entries)
            line_printf('%s\n    covers: %s\n', entries(i).ref, entries(i).covers);
        end
    end
    clear entries
end
end
