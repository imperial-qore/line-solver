function [tex, sys] = exportODEs(self, filename, notation)
% [TEX, SYS] = EXPORTODES(FILENAME, NOTATION)
% Export the system of ODEs integrated by the mean-field methods of
% SolverFLD as a standalone LaTeX document, in a symbolic form that is
% both human and machine readable. The exported system corresponds to
% the solver method set in the solver options (default/matrix, pnorm,
% closing, statedep, softmin).
%
% Input:
%  filename - path of the .tex file to write; if empty or omitted the
%             LaTeX source is returned without writing a file
%  notation - 'scalar' (default): one expanded ODE per state variable
%             'matrix': compact matrix notation, dx/dt = W'*theta(x) +
%             lambda for the matrix/pnorm methods and dx/dt = J*r(x) for
%             the closing/statedep/softmin methods
%
% Output:
%  tex - LaTeX source (char)
%  sys - structural description of the ODE system (see solver_fluid_symodes)
%
% The document header also carries the state-variable mapping as LaTeX
% comments (lines starting with % STATE) for machine parsing.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 2
    filename = '';
end
if nargin < 3
    notation = 'scalar';
end
switch notation
    case {'scalar','matrix'}
        % valid
    otherwise
        line_error(mfilename, sprintf('Unknown notation ''%s''. Valid notations: scalar, matrix.', notation));
end

% lang='cpp' exports the drift from line-cli (-a odes), with --notation
% reaching the C++ exporter so 'matrix' returns the matrix document and not
% the scalar one under another name. SYS is the structural description
% solver_fluid_symodes builds on this side; the CLI sends the document only,
% so SYS is empty rather than reconstructed from another engine's numbers.
if isfield(self.options,'lang') && strcmp(self.options.lang,'cpp')
    tex = CPPLINE.exportODEs(self.name, self.model, self.options, notation);
    sys = struct([]);
    if ~isempty(filename)
        fid = fopen(filename, 'w');
        if fid < 0
            line_error(mfilename, sprintf('could not open ''%s'' for writing.', filename));
        end
        fprintf(fid, '%s', tex);
        fclose(fid);
    end
    return
end

options = self.getOptions;
sn = self.getStruct;
% SSA draws a sample path and the fluid ODEs read mu*phi as a flow, so their
% surrogate must be a genuine phase-type: a matrix exponential has neither.
options.config.phfit = 'ph';
sn = sn_nonmarkov_toph(sn, options);
sys = solver_fluid_symodes(sn, options);
sys.x0 = build_x0(sys, sn, options);

modelName = self.model.getName;
n = sys.nstates;

%% header comments (machine readable)
L = {};
L{end+1} = '% Mean-field fluid ODE system exported by LINE SolverFLD';
L{end+1} = sprintf('%% model: %s', modelName);
L{end+1} = sprintf('%% method: %s', sys.method);
switch sys.form
    case 'W'
        L{end+1} = '% form: dx/dt = W^T*theta(x) + lambda';
    case 'J'
        L{end+1} = '% form: dx/dt = J*r(x)';
end
L{end+1} = sprintf('%% notation: %s', notation);
L{end+1} = sprintf('%% nstates: %d', n);
if strcmp(sys.form,'J')
    L{end+1} = sprintf('%% nevents: %d', sys.nevents);
end
for s = 1:n
    L{end+1} = sprintf('%% STATE %d station=%s class=%s phase=%d', s, ...
        sys.stationNames{sys.stateStation(s)}, sys.classNames{sys.stateClass(s)}, sys.statePhase(s)); %#ok<AGROW>
end
if strcmp(sys.form,'J')
    for e = 1:sys.nevents
        L{end+1} = sprintf('%% EVENT %d var=%d type=%s coeff=%.15g', e, ...
            sys.eventVar(e), sys.factorType{e}, sys.coeff(e)); %#ok<AGROW>
    end
end

%% preamble
L{end+1} = '\documentclass{article}';
L{end+1} = '\usepackage{amsmath}';
L{end+1} = '\usepackage[margin=2.5cm]{geometry}';
L{end+1} = '\allowdisplaybreaks';
L{end+1} = '\setcounter{MaxMatrixCols}{500}';
L{end+1} = '\begin{document}';
L{end+1} = '\section*{Mean-field fluid ODE system}';
L{end+1} = sprintf('\\noindent Model: \\texttt{%s}. Solver: \\texttt{SolverFLD}, method \\texttt{%s}, %s notation.', ...
    texesc(modelName), texesc(sys.method), notation);
switch sys.form
    case 'W'
        L{end+1} = sprintf('The system has %d state variables and reads $\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = W^{\\top}\\theta(\\mathbf{x}) + \\boldsymbol{\\lambda}$.', n);
    case 'J'
        L{end+1} = sprintf('The system has %d state variables and %d events and reads $\\frac{\\mathrm{d}\\mathbf{x}}{\\mathrm{d}t} = J\\,r(\\mathbf{x})$.', n, sys.nevents);
end

%% state variable legend
L{end+1} = '\subsection*{State variables}';
L{end+1} = 'Each state variable $x_{s}$ is the mean number of jobs of a class in a service phase at a station:';
L{end+1} = '\begin{center}';
chunk = 48;
for s0 = 1:chunk:n
    s1 = min(n, s0+chunk-1);
    L{end+1} = '\begin{tabular}{rlll}'; %#ok<AGROW>
    L{end+1} = '\hline'; %#ok<AGROW>
    L{end+1} = '$s$ & station & class & phase\\'; %#ok<AGROW>
    L{end+1} = '\hline'; %#ok<AGROW>
    for s = s0:s1
        L{end+1} = sprintf('%d & \\texttt{%s} & \\texttt{%s} & %d\\\\', s, ...
            texesc(sys.stationNames{sys.stateStation(s)}), ...
            texesc(sys.classNames{sys.stateClass(s)}), sys.statePhase(s)); %#ok<AGROW>
    end
    L{end+1} = '\hline'; %#ok<AGROW>
    L{end+1} = '\end{tabular}'; %#ok<AGROW>
    if s1 < n
        L{end+1} = '\par\medskip'; %#ok<AGROW>
    end
end
L{end+1} = '\end{center}';

% station legend
usedStations = unique(sys.stateStation(:))';
L{end+1} = '\begin{center}';
L{end+1} = '\begin{tabular}{rlll}';
L{end+1} = '\hline';
L{end+1} = '$i$ & station & scheduling & $S_{i}$\\';
L{end+1} = '\hline';
for i = usedStations
    L{end+1} = sprintf('%d & \\texttt{%s} & %s & $%s$\\\\', i, ...
        texesc(sys.stationNames{i}), texesc(sys.schedNames{i}), fmtnum(sys.S(i))); %#ok<AGROW>
end
L{end+1} = '\hline';
L{end+1} = '\end{tabular}';
L{end+1} = '\end{center}';

%% definitions of station-level auxiliary functions
[T, varFactor, constTerm] = build_terms(sys);
defs = build_defs(sys, varFactor);
if ~isempty(defs)
    L{end+1} = '\subsection*{Definitions}';
    L{end+1} = '\begin{align*}';
    for d = 1:length(defs)
        L{end+1} = defs{d}; %#ok<AGROW>
    end
    L{end+1} = '\end{align*}';
end

%% ODE system
switch notation
    case 'scalar'
        L{end+1} = '\subsection*{ODE system (scalar notation)}';
        L{end+1} = '\begin{align}';
        for s = 1:n
            L{end+1} = render_equation(sys, s, T, varFactor, constTerm, s==n); %#ok<AGROW>
        end
        L{end+1} = '\end{align}';
    case 'matrix'
        L{end+1} = '\subsection*{ODE system (matrix notation)}';
        switch sys.form
            case 'W'
                haveLambda = any(sys.Alambda ~= 0);
                if haveLambda
                    L{end+1} = '\begin{equation}';
                    L{end+1} = '\frac{\mathrm{d}\mathbf{x}}{\mathrm{d}t} = W^{\top}\,\theta(\mathbf{x}) + \boldsymbol{\lambda}';
                    L{end+1} = '\end{equation}';
                else
                    L{end+1} = '\begin{equation}';
                    L{end+1} = '\frac{\mathrm{d}\mathbf{x}}{\mathrm{d}t} = W^{\top}\,\theta(\mathbf{x})';
                    L{end+1} = '\end{equation}';
                end
                L{end+1} = 'with $\theta_{s}(\mathbf{x})$ given componentwise by';
                L{end+1} = '\begin{equation*}';
                L{end+1} = ['\theta(\mathbf{x}) = ', render_theta_vector(sys, varFactor)];
                L{end+1} = '\end{equation*}';
                L{end+1} = 'and';
                L{end+1} = '\begin{equation*}';
                L{end+1} = ['W^{\top} = ', render_num_matrix(sys.W')];
                L{end+1} = '\end{equation*}';
                if haveLambda
                    L{end+1} = '\begin{equation*}';
                    L{end+1} = ['\boldsymbol{\lambda} = ', render_num_vector(sys.Alambda), '^{\top}'];
                    L{end+1} = '\end{equation*}';
                end
            case 'J'
                L{end+1} = '\begin{equation}';
                L{end+1} = '\frac{\mathrm{d}\mathbf{x}}{\mathrm{d}t} = J\,r(\mathbf{x})';
                L{end+1} = '\end{equation}';
                L{end+1} = 'with stoichiometry matrix';
                L{end+1} = '\begin{equation*}';
                L{end+1} = ['J = ', render_num_matrix(sys.J)];
                L{end+1} = '\end{equation*}';
                L{end+1} = 'and event rate functions';
                L{end+1} = '\begin{align*}';
                for e = 1:sys.nevents
                    fd = sys.factorData{e};
                    fstr = factor_tex(sys.eventVar(e), sys.factorType{e}, fd);
                    L{end+1} = sprintf('r_{%d}(\\mathbf{x}) &= %s%s', e, ...
                        term_tex(sys.coeff(e), fstr), tern(e<sys.nevents,'\\','')); %#ok<AGROW>
                end
                L{end+1} = '\end{align*}';
        end
end

%% initial condition
if ~isempty(sys.x0)
    L{end+1} = '\subsection*{Initial condition}';
    L{end+1} = '\begin{equation*}';
    L{end+1} = ['\mathbf{x}(0) = ', render_num_vector(sys.x0), '^{\top}'];
    L{end+1} = '\end{equation*}';
end

%% remarks
L{end+1} = '\subsection*{Remarks}';
L{end+1} = '\begin{itemize}';
L{end+1} = sprintf('\\item For each station $i$, $n_{i}(\\mathbf{x})$ denotes the total mass at the station and $S_{i}$ the number of servers (infinite-server stations use the closed job population, $%s$ denotes infinity).', '\infty');
L{end+1} = '\item The numerical solver regularizes vanishing denominators with a small positive constant; these regularizations are omitted here.';
if strcmp(sys.form,'J') && any(strcmp('fcfsw', sys.factorType) | strcmp('fcfsws', sys.factorType))
    L{end+1} = '\item At FCFS stations, the mean phase residence times $w_{u} = -1/[D_{0}]_{kk}$ weight the backlog $\hat{n}_{i}$; the factors $w_{u}$ of the departing phases are folded into the rate coefficients.';
end
if strcmp(sys.form,'J') && any(strcmp('dpsmin', sys.factorType))
    L{end+1} = '\item At DPS stations, weights are normalized to sum to one and the weight $w_{ir}$ of the departing class is folded into the rate coefficient; the class shares $w_{ir}x/\tilde{n}_{i}$ divide the station capacity $\min(n_{i},S_{i})$, so they sum to one whenever the station is busy.';
end
if any(sys.sched(unique(sys.stateStation)) == SchedStrategy.FCFS) && any(strcmp(sys.method, {'matrix','closing'}))
    L{end+1} = '\item For FCFS stations with non-exponential service, the solver may iteratively re-fit the service distributions (non-exponential approximation); the exported system uses the nominal model parameters.';
end
if isfield(options.config,'hide_immediate') && options.config.hide_immediate
    L{end+1} = '\item \texttt{hide\_immediate} is enabled in the solver options: the numerical integration may further eliminate immediate transitions by state-space reduction; the exported system is the unreduced one.';
end
L{end+1} = '\end{itemize}';
L{end+1} = '\end{document}';

tex = strjoin(L, newline);
tex = sprintf('%s\n', tex);

if ~isempty(filename)
    fid = fopen(filename, 'w');
    if fid == -1
        line_error(mfilename, sprintf('Cannot open file ''%s'' for writing.', filename));
    end
    fprintf(fid, '%s', tex);
    fclose(fid);
end
end

%% helpers

function x0 = build_x0(sys, sn, options)
% initial condition in the exported state space
init_sol = solver_fluid_initsol(sn, options);
init_sol = init_sol(:);
switch sys.form
    case 'J'
        x0 = init_sol;
    case 'W'
        M = sn.nstations;
        K = sn.nclasses;
        nphases = sn.phases;
        x0_build = [];
        state = 0;
        init_idx = 0;
        for ist = 1:M
            for r = 1:K
                if nphases(ist,r) == 0
                    state = state + 1;
                    x0_build(state,1) = 0; %#ok<AGROW>
                else
                    for k = 1:nphases(ist,r)
                        state = state + 1;
                        if isnan(sn.rates(ist,r))
                            x0_build(state,1) = 0; %#ok<AGROW>
                        else
                            init_idx = init_idx + 1;
                            x0_build(state,1) = init_sol(init_idx); %#ok<AGROW>
                        end
                    end
                end
            end
        end
        x0 = x0_build(sys.keep);
        x0(sys.isSource) = 0;
end
end

function [T, varFactor, constTerm] = build_terms(sys)
% T(s,v): constant coefficient of the term driven by state variable v in
% the equation of state s; varFactor{v}: factor descriptor of variable v;
% constTerm(s): additive constant.
n = sys.nstates;
T = zeros(n,n);
varFactor = cell(n,1);
switch sys.form
    case 'W'
        T = sys.W';
        T(:, sys.isSource) = 0; % theta of Source states is identically zero
        for v = 1:n
            if ~sys.isSource(v)
                i = sys.stateStation(v);
                % sys.S holds the population at an INF station, never Inf
                if sys.isInfStation(i)
                    varFactor{v} = struct('type','lin','station',i,'class',sys.stateClass(v));
                else
                    varFactor{v} = struct('type',sys.smoothing,'station',i,'class',sys.stateClass(v));
                end
            end
        end
        constTerm = sys.Alambda;
    case 'J'
        for e = 1:sys.nevents
            v = sys.eventVar(e);
            T(:,v) = T(:,v) + sys.J(:,e) * sys.coeff(e);
            if isempty(varFactor{v})
                fd = sys.factorData{e};
                fd.type = sys.factorType{e};
                varFactor{v} = fd;
            end
        end
        constTerm = zeros(n,1);
end
end

function defs = build_defs(sys, varFactor)
% station-level auxiliary definitions used by the factors
defs = {};
n = sys.nstates;
M = length(sys.stationNames);
needN = false(M,1);
needNT = false(M,1); % ntilde (DPS-weighted mass)
needNH = false(M,1); % nhat (FCFS backlog)
gdef = cell(M,1);
for v = 1:n
    f = varFactor{v};
    if isempty(f)
        continue
    end
    i = f.station;
    switch f.type
        case 'lin'
            % no auxiliary quantity
        case 'min'
            needN(i) = true;
            gdef{i} = sprintf('g_{%d}(\\mathbf{x}) &= \\frac{\\min(n_{%d}(\\mathbf{x}),\\, %s)}{n_{%d}(\\mathbf{x})}', i, i, fmtnum(sys.S(i)), i);
        case 'pnorm'
            needN(i) = true;
            gdef{i} = sprintf('g_{%d}(\\mathbf{x}) &= \\Bigl(1 + \\bigl(n_{%d}(\\mathbf{x})/%s\\bigr)^{%s}\\Bigr)^{-1/%s}', ...
                i, i, fmtnum(sys.S(i)), fmtnum(sys.pstar(i)), fmtnum(sys.pstar(i)));
        case 'dpsmin'
            needN(i) = true;
            needNT(i) = true;
            gdef{i} = sprintf('g_{%d}(\\mathbf{x}) &= \\frac{\\min(n_{%d}(\\mathbf{x}),\\, %s)}{\\tilde{n}_{%d}(\\mathbf{x})}', i, i, fmtnum(sys.S(i)), i);
        case 'dpspw'
            needN(i) = true;
            needNT(i) = true;
            if isempty(gdef{i})
                gdef{i} = '';
            end
            % one definition per enabled class at the station, built below
        case {'fcfsw','fcfsws'}
            needN(i) = true;
            needNH(i) = true;
            if strcmp(f.type,'fcfsw')
                gdef{i} = sprintf('g_{%d}(\\mathbf{x}) &= \\frac{\\min(n_{%d}(\\mathbf{x}),\\, %s)}{\\hat{n}_{%d}(\\mathbf{x})}', i, i, fmtnum(sys.S(i)), i);
            else
                gdef{i} = sprintf('g_{%d}(\\mathbf{x}) &= \\frac{\\mathrm{softmin}\\bigl(n_{%d}(\\mathbf{x}),\\, %s\\bigr)}{\\hat{n}_{%d}(\\mathbf{x})}', i, i, fmtnum(sys.S(i)), i);
            end
    end
end
% n_i definitions
for i = 1:M
    if needN(i)
        vlist = find(sys.stateStation == i);
        defs{end+1} = sprintf('n_{%d}(\\mathbf{x}) &= %s\\\\', i, strjoin(arrayfun(@(v) sprintf('x_{%d}',v), vlist(:)', 'UniformOutput', false), ' + ')); %#ok<AGROW>
    end
end
% ntilde_i definitions (DPS-weighted station mass)
for i = 1:M
    if needNT(i)
        parts = {};
        K = length(sys.classNames);
        for r = 1:K
            vlist = find(sys.stateStation == i & sys.stateClass == r);
            if ~isempty(vlist)
                parts{end+1} = sprintf('%s\\,(%s)', fmtnum(sys.dpsw(i,r)), strjoin(arrayfun(@(v) sprintf('x_{%d}',v), vlist(:)', 'UniformOutput', false), ' + ')); %#ok<AGROW>
            end
        end
        defs{end+1} = sprintf('\\tilde{n}_{%d}(\\mathbf{x}) &= %s\\\\', i, strjoin(parts, ' + ')); %#ok<AGROW>
    end
end
% nhat_i definitions (FCFS backlog weighted by mean phase residence times)
for i = 1:M
    if needNH(i)
        vlist = find(sys.stateStation == i);
        parts = arrayfun(@(v) sprintf('%s\\,x_{%d}', fmtnum(sys.fcfsPhaseW(v)), v), vlist(:)', 'UniformOutput', false);
        defs{end+1} = sprintf('\\hat{n}_{%d}(\\mathbf{x}) &= %s\\\\', i, strjoin(parts, ' + ')); %#ok<AGROW>
    end
end
% g definitions
for i = 1:M
    if ~isempty(gdef{i})
        defs{end+1} = sprintf('%s\\\\', gdef{i}); %#ok<AGROW>
    end
end
% per-class piecewise DPS definitions
for v = 1:n
    f = varFactor{v};
    if ~isempty(f) && strcmp(f.type,'dpspw')
        i = f.station;
        r = f.class;
        d = sprintf('g_{%d,%d}(\\mathbf{x}) &= \\begin{cases} 1 & n_{%d}(\\mathbf{x}) \\le %s\\\\ \\dfrac{%s}{\\tilde{n}_{%d}(\\mathbf{x})} & n_{%d}(\\mathbf{x}) > %s \\end{cases}\\\\', ...
            i, r, i, fmtnum(sys.S(i)), fmtnum(sys.S(i)*sys.dpsw(i,r)), i, i, fmtnum(sys.S(i)));
        if ~any(strcmp(defs, d))
            defs{end+1} = d; %#ok<AGROW>
        end
    end
end
% softmin definition
if isfield(sys,'alpha') && any(cellfun(@(f) ~isempty(f) && strcmp(f.type,'fcfsws'), varFactor))
    defs{end+1} = sprintf('\\mathrm{softmin}(a,b) &= \\frac{a\\,e^{-\\alpha a} + b\\,e^{-\\alpha b}}{e^{-\\alpha a} + e^{-\\alpha b}}, \\qquad \\alpha = %s\\\\', fmtnum(sys.alpha));
end
% strip the trailing line break of the last definition
if ~isempty(defs)
    defs{end} = regexprep(defs{end}, '\\\\$', '');
end
end

function s = render_equation(sys, sidx, T, varFactor, constTerm, isLast)
% one align row for the ODE of state variable sidx
terms = {};
for v = 1:sys.nstates
    c = T(sidx,v);
    if c ~= 0
        fstr = factor_tex(v, varFactor{v}.type, varFactor{v});
        terms{end+1} = {c, fstr}; %#ok<AGROW>
    end
end
if constTerm(sidx) ~= 0
    terms{end+1} = {constTerm(sidx), ''};
end
if isempty(terms)
    rhs = '0';
else
    parts = {};
    for k = 1:length(terms)
        c = terms{k}{1};
        body = term_tex(abs(c), terms{k}{2});
        if k == 1
            if c < 0
                parts{end+1} = ['-', body]; %#ok<AGROW>
            else
                parts{end+1} = body; %#ok<AGROW>
            end
        else
            if c < 0
                parts{end+1} = [' - ', body]; %#ok<AGROW>
            else
                parts{end+1} = [' + ', body]; %#ok<AGROW>
            end
        end
        % break long equations every four terms
        if mod(k,4) == 0 && k < length(terms)
            parts{end+1} = sprintf('\\nonumber\\\\\n&\\quad '); %#ok<AGROW>
        end
    end
    rhs = strjoin(parts, '');
end
s = sprintf('\\frac{\\mathrm{d}x_{%d}}{\\mathrm{d}t} &= %s%s', sidx, rhs, tern(~isLast,'\\',''));
end

function fstr = factor_tex(v, ftype, fdata)
% LaTeX of the state-dependent factor of a term driven by variable v
switch ftype
    case 'lin'
        fstr = sprintf('x_{%d}', v);
    case {'min','pnorm','dpsmin','fcfsw','fcfsws'}
        fstr = sprintf('x_{%d}\\,g_{%d}(\\mathbf{x})', v, fdata.station);
    case 'dpspw'
        fstr = sprintf('x_{%d}\\,g_{%d,%d}(\\mathbf{x})', v, fdata.station, fdata.class);
    case 'ext1'
        if isempty(fdata.others)
            fstr = ''; % single-phase source class: constant unit mass
        else
            fstr = sprintf('\\bigl(1 - %s\\bigr)', strjoin(arrayfun(@(u) sprintf('x_{%d}',u), fdata.others(:)', 'UniformOutput', false), ' - '));
        end
end
end

function s = term_tex(c, fstr)
% LaTeX of a term with positive coefficient c and factor fstr
if isempty(fstr)
    s = fmtnum(c);
elseif c == 1
    s = fstr;
else
    s = sprintf('%s\\,%s', fmtnum(c), fstr);
end
end

function s = render_theta_vector(sys, varFactor)
% componentwise theta vector for the matrix notation of the W form
rows = cell(sys.nstates,1);
for v = 1:sys.nstates
    if sys.isSource(v)
        rows{v} = '0';
    else
        rows{v} = factor_tex(v, varFactor{v}.type, varFactor{v});
    end
end
s = sprintf('\\begin{bmatrix} %s \\end{bmatrix}', strjoin(rows, ' \\\\ '));
end

function s = render_num_matrix(A)
% numeric matrix as a LaTeX bmatrix
[m, ncol] = size(A); %#ok<ASGLU>
rows = cell(m,1);
for r = 1:m
    rows{r} = strjoin(arrayfun(@fmtnum, A(r,:), 'UniformOutput', false), ' & ');
end
body = strjoin(rows, ' \\\\ ');
if max(size(A)) > 12
    s = sprintf('{\\scriptsize\\begin{bmatrix} %s \\end{bmatrix}}', body);
else
    s = sprintf('\\begin{bmatrix} %s \\end{bmatrix}', body);
end
end

function s = render_num_vector(v)
% numeric column vector as a transposed LaTeX row
s = sprintf('\\begin{pmatrix} %s \\end{pmatrix}', strjoin(arrayfun(@fmtnum, v(:)', 'UniformOutput', false), ' & '));
end

function s = fmtnum(v)
% compact LaTeX-safe number formatting
if isinf(v)
    if v > 0
        s = '\infty';
    else
        s = '-\infty';
    end
elseif v == round(v) && abs(v) < 1e15
    s = sprintf('%d', v);
else
    s = sprintf('%.8g', v);
end
end

function s = texesc(s)
% escape LaTeX special characters in identifiers
s = regexprep(s, '([_%&#])', '\\$1');
end

function out = tern(cond, a, b)
if cond
    out = a;
else
    out = b;
end
end
