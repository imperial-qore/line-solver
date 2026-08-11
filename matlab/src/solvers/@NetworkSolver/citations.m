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

% the solver family, both as a token of its own (for the paradigms whose
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
end

% the percentile method of the last getPerctRespT call, if any
if ~isempty(self.lastPerctMethod)
    tokens{end+1} = self.lastPerctMethod;
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
