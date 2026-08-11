function [SensTable, sens] = getSensitivityTable(self, varargin)
% GETSENSITIVITYTABLE Layer-wise performance sensitivities of a layered network.
%
% [SENSTABLE, SENS] = GETSENSITIVITYTABLE(SELF) solves the layered model and
% then delegates to each layer solver, returning the concatenation of the layer
% tables with a leading Layer column. Every row is therefore a (Layer, Station,
% JobClass) triple carrying the derivative of that row's mean measures with
% respect to that station-class service RATE in that layer:
%   dTput_dRate, dRespT_dRate, dQLen_dRate, dUtil_dRate.
%
% [...] = GETSENSITIVITYTABLE(SELF, 'method', M, 'step', H, 'scheme', S) passes
% the options through to the layer solvers unchanged, with the same meaning as
% in NetworkSolver/getSensitivityTable: each layer independently takes the
% analytic branch where its own solver supports it and the model is in scope,
% and finite differences otherwise. SENS is a cell array of the layer second
% outputs, indexed by layer.
%
% IMPORTANT, on what these derivatives mean. Each entry is a derivative WITHIN
% ITS LAYER, taken with the layer parameters that the fixed point produced held
% fixed. It is a partial derivative of the layer submodel, not the total
% derivative of the layered model: perturbing a host demand in one layer moves
% the think times, populations and service rates of the other layers through
% the fixed-point map, and that indirect term is not included here. The layer
% table is the right object for attributing a bottleneck inside a layer, and
% the wrong one for predicting the effect of a parameter change on the solved
% layered model. For the latter, finite-difference the LayeredNetwork itself.

if (GlobalConstants.DummyMode)
    SensTable = [];
    sens = {};
    return
end

% see _kb/06-solver-catalog.md (LN section) for rationale
if ~isempty(self.obj)
    [SensTable, sens] = i_javaSensitivityTable(self, varargin{:});
    return
elseif ~isempty(self.pyMode) && self.pyMode
    [SensTable, sens] = PYLINE.getLNSensitivityTable(self.model, self.options, varargin{:});
    return
end

% see _kb/06-solver-catalog.md (LN section) for rationale
if isempty(self.results)
    self.iterate();
end

E = self.getNumberOfModels();
Layer = {}; Station = {}; JobClass = {};
dTput = []; dRespT = []; dQLen = []; dUtil = [];
sens = cell(1, E);
methods = cell(1, E);

for e = 1:E
    solver = self.solvers{e};
    if isempty(solver)
        continue
    end
    [T, sens{e}] = solver.getSensitivityTable(varargin{:});
    methods{e} = T.Properties.UserData.method;
    layerName = self.ensemble{e}.getName();
    for r = 1:height(T)
        Layer{end+1, 1}    = layerName;      %#ok<AGROW>
        Station{end+1, 1}  = T.Station{r};   %#ok<AGROW>
        JobClass{end+1, 1} = T.JobClass{r};  %#ok<AGROW>
        dTput(end+1, 1)  = T.dTput_dRate(r);  %#ok<AGROW>
        dRespT(end+1, 1) = T.dRespT_dRate(r); %#ok<AGROW>
        dQLen(end+1, 1)  = T.dQLen_dRate(r);  %#ok<AGROW>
        dUtil(end+1, 1)  = T.dUtil_dRate(r);  %#ok<AGROW>
    end
end

SensTable = table(Layer, Station, JobClass, dTput, dRespT, dQLen, dUtil, ...
    'VariableNames', {'Layer', 'Station', 'JobClass', 'dTput_dRate', ...
    'dRespT_dRate', 'dQLen_dRate', 'dUtil_dRate'});

% One branch label per layer, plus a summary that is 'mixed' when the layers
% did not all take the same branch.
present = methods(~cellfun(@isempty, methods));
if isempty(present)
    summary = '';
elseif all(strcmp(present, present{1}))
    summary = present{1};
else
    summary = 'mixed';
end
SensTable.Properties.UserData = struct('method', summary, 'layerMethods', {methods});
end

function [SensTable, sens] = i_javaSensitivityTable(self, varargin)
% I_JAVASENSITIVITYTABLE  lang='java': delegate to the JAR SolverLN, which owns
% the layer ensemble in this mode, and marshal its LayeredNetworkSensitivityTable
% back into the same MATLAB table and UserData that the native branch returns.
opt = i_parseSensOptions(varargin{:});
step = opt.step;
if isempty(step)
    step = NaN; % NaN selects the JAR-side default, as in the Java signature
end
jt = self.obj.getSensitivityTable(opt.method, step, opt.scheme);

Layer    = i_javaStringList(jt.getLayerNames());
Station  = i_javaStringList(jt.getStationNames());
JobClass = i_javaStringList(jt.getClassNames());
dTput  = i_javaDoubleList(jt.getDTput());
dRespT = i_javaDoubleList(jt.getDRespT());
dQLen  = i_javaDoubleList(jt.getDQLen());
dUtil  = i_javaDoubleList(jt.getDUtil());

SensTable = table(Layer, Station, JobClass, dTput, dRespT, dQLen, dUtil, ...
    'VariableNames', {'Layer', 'Station', 'JobClass', 'dTput_dRate', ...
    'dRespT_dRate', 'dQLen_dRate', 'dUtil_dRate'});

methods = i_javaStringList(jt.getLayerMethods())';
jsens = jt.getLayerSens();
sens = cell(1, numel(methods));
for e = 1:numel(sens)
    sens{e} = i_javaSensStruct(jsens.get(e-1));
end
SensTable.Properties.UserData = struct('method', char(jt.getMethod()), ...
    'layerMethods', {methods});
end

function opt = i_parseSensOptions(varargin)
% I_PARSESENSOPTIONS  Same name-value contract as
% @NetworkSolver/getSensitivityTable, parsed here because the bridged branches
% need the values rather than passing varargin straight through.
opt = struct('method', 'auto', 'step', [], 'scheme', 'forward');
if mod(numel(varargin), 2) ~= 0
    line_error(mfilename, 'Options must be given as name-value pairs.');
end
for a = 1:2:numel(varargin)
    name = lower(varargin{a});
    if ~isfield(opt, name)
        line_error(mfilename, sprintf('Unknown option ''%s''.', varargin{a}));
    end
    opt.(name) = varargin{a+1};
end
if ischar(opt.method), opt.method = lower(opt.method); end
if ischar(opt.scheme), opt.scheme = lower(opt.scheme); end
if ~ismember(opt.method, {'auto', 'exact', 'fd'})
    line_error(mfilename, 'The method must be one of ''auto'', ''exact'', ''fd''.');
end
if ~ismember(opt.scheme, {'forward', 'central'})
    line_error(mfilename, 'The scheme must be ''forward'' or ''central''.');
end
end

function c = i_javaStringList(jlist)
% I_JAVASTRINGLIST  java.util.List<String> -> column cell of char. A null entry
% (a layer with no solver) becomes '' so that the cellfun(@isempty) tests and
% the strcmp summary above behave as on the native branch.
n = jlist.size();
c = cell(n, 1);
for i = 1:n
    e = jlist.get(i-1);
    if isempty(e)
        c{i} = '';
    else
        c{i} = char(e);
    end
end
end

function v = i_javaDoubleList(jlist)
% I_JAVADOUBLELIST  java.util.List<Double> -> column double vector.
n = jlist.size();
v = zeros(n, 1);
for i = 1:n
    v(i) = double(jlist.get(i-1));
end
end

function s = i_javaSensStruct(jsens)
% I_JAVASENSSTRUCT  jline.io.Ret.pfqnSens -> the pfqn_sens struct that the
% native MATLAB branch returns as its second output. Empty for a layer that
% took the finite-difference branch, which carries no analytic Jacobian.
if isempty(jsens)
    s = [];
    return
end
s = struct();
s.X = JLINE.from_jline_matrix(jsens.X);
s.Q = JLINE.from_jline_matrix(jsens.Q);
s.U = JLINE.from_jline_matrix(jsens.U);
s.R = JLINE.from_jline_matrix(jsens.R);
s.dX = JLINE.from_jline_matrix(jsens.dX);
s.dQ = i_javaMatrixArray3(jsens.dQ);
s.dU = i_javaMatrixArray3(jsens.dU);
s.dR = i_javaMatrixArray3(jsens.dR);
% The JAR keeps the parameter descriptors as three parallel int arrays with
% type 0 = L / 1 = Z and station -1 for Z; pfqn_sens.m returns one 1 x P struct
% array with .type 'L'/'Z', .station (0 for Z) and .class. Convert, so that a
% caller reading sens{e}.params(p) does not have to know which backend ran.
jType = double(jsens.paramType(:))';
jStation = double(jsens.paramStation(:))';
jClass = double(jsens.paramClass(:))';
P = numel(jType);
params = struct('type', cell(1, P), 'station', cell(1, P), 'class', cell(1, P));
for p = 1:P
    if jType(p) == 0
        params(p).type = 'L';
        params(p).station = jStation(p) + 1;   % 0-based station -> 1-based
    else
        params(p).type = 'Z';
        params(p).station = 0;
    end
    params(p).class = jClass(p) + 1;           % 0-based class -> 1-based
end
s.params = params;
if ~isempty(jsens.QVar)
    s.QVar = JLINE.from_jline_matrix(jsens.QVar);
end
if ~isempty(jsens.QTotVar)
    s.QTotVar = JLINE.from_jline_matrix(jsens.QTotVar);
end
s.QCovAsym = double(jsens.QCovAsym);
if ~isempty(jsens.QCov)
    % QCov[i][r] is an M x R matrix of Cov[n(i,r),n(j,s)]; pfqn_sens.m returns
    % the same content as one M x R x M x R array, so index it that way here.
    M = numel(jsens.QCov);
    firstRow = jsens.QCov(1);
    Rn = size(JLINE.from_jline_matrix(firstRow(1)), 2);
    s.QCov = zeros(M, Rn, M, Rn);
    for i = 1:M
        row = jsens.QCov(i);
        for r = 1:Rn
            s.QCov(i, r, :, :) = JLINE.from_jline_matrix(row(r));
        end
    end
end
end

function A = i_javaMatrixArray3(jarr)
% I_JAVAMATRIXARRAY3  Matrix[] of P entries, each M x R -> M x R x P array,
% matching the shape pfqn_sens returns in MATLAB (dQ(ist,c,p)).
P = numel(jarr);
if P == 0
    A = [];
    return
end
first = JLINE.from_jline_matrix(jarr(1));
A = zeros(size(first, 1), size(first, 2), P);
A(:, :, 1) = first;
for p = 2:P
    A(:, :, p) = JLINE.from_jline_matrix(jarr(p));
end
end
