function [ok, msg] = ctmc_memory_gate(logNstates, force, verbose, safetyFraction)
% CTMC_MEMORY_GATE  Hardware-aware, profiling-calibrated CTMC memory pre-gate.
%
% [ok,msg] = ctmc_memory_gate(logNstates, force, verbose, safetyFraction)
%
% Decides whether a CTMC steady-state solve of a state space of worst-case
% size exp(logNstates) is safe on the current host. The budget is a fraction
% of the memory actually available (see lineGetAvailableMemory); the per-state
% cost is calibrated by profiling the local sparse LU once and caching the
% fitted power law per machine in tempdir.
%
% ok is false only when the predicted peak footprint exceeds the budget and
% force is not set. This replaces the historical hard-coded threshold.

if nargin < 2 || isempty(force); force = false; end
if nargin < 3 || isempty(verbose); verbose = false; end
if nargin < 4 || isempty(safetyFraction); safetyFraction = 0.6; end

avail  = lineGetAvailableMemory();
budget = safetyFraction * avail;
calib  = ctmc_get_calibration(verbose);

logPred   = log(calib.alpha_mem) + calib.beta_mem * logNstates;
logBudget = log(max(budget, 1));
predGB    = exp(min(logPred, 700)) / 1024^3;
budgetGB  = budget / 1024^3;

ok = true;
msg = '';
if logPred > logBudget
    msg = sprintf(['CTMC predicted peak memory ~%.2f GB exceeds the safe ' ...
        'budget ~%.2f GB (%.0f%% of %.2f GB available). Reduce the state ' ...
        'space (e.g. lower ''cutoff''), use another solver (MVA/NC/FLD), or ' ...
        'set force=true to override.'], predGB, budgetGB, ...
        100*safetyFraction, avail/1024^3);
    if ~force
        ok = false;
        return
    end
    if verbose
        line_printf('Warning (forced): %s\n', msg);
    end
elseif verbose && logPred > log(max(0.5*budget, 1))
    line_printf('CTMC predicted peak memory ~%.2f GB (budget ~%.2f GB).\n', ...
        predGB, budgetGB);
end
end

% ------------------------------------------------------------------------
function calib = ctmc_get_calibration(verbose)
% Return calibrated power-law coefficients {alpha_mem,beta_mem,alpha_t,beta_t}
% for factor bytes and factorization time, using a per-machine tempdir cache.
BYTES_PER_NZ = 16;
G1 = 40; G2 = 80;
FALLBACK_ALPHA = BYTES_PER_NZ * 8;
FALLBACK_BETA  = 1.3;

sig = ctmc_machine_signature();
cachefile = fullfile(tempdir, 'line_ctmc_calib_matlab.json');

if exist(cachefile, 'file')
    try
        txt = fileread(cachefile);
        data = jsondecode(txt);
        if isfield(data, 'sig') && strcmp(data.sig, sig)
            calib = data;
            return
        end
    catch
    end
end

try
    [n1, b1, t1] = ctmc_profile_point(G1, BYTES_PER_NZ);
    [n2, b2, t2] = ctmc_profile_point(G2, BYTES_PER_NZ);
    [am, bm] = ctmc_fit_power_law(n1, b1, n2, b2);
    [at, bt] = ctmc_fit_power_law(n1, max(t1,1e-9), n2, max(t2,1e-9));
    calib = struct('sig', sig, 'alpha_mem', am, 'beta_mem', bm, ...
        'alpha_t', at, 'beta_t', bt, 'timestamp', posixtime(datetime('now')));
    try
        fid = fopen(cachefile, 'w');
        if fid > 0
            fwrite(fid, jsonencode(calib));
            fclose(fid);
        end
    catch
    end
    if verbose
        line_printf('CTMC calibration: bytes ~ %.3g*N^%.3f\n', am, bm);
    end
catch ME
    if verbose
        line_printf('CTMC calibration failed (%s); using fallback model\n', ME.message);
    end
    calib = struct('sig', sig, 'alpha_mem', FALLBACK_ALPHA, ...
        'beta_mem', FALLBACK_BETA, 'alpha_t', 0, 'beta_t', 1, 'timestamp', posixtime(datetime('now')));
end
end

% ------------------------------------------------------------------------
function sig = ctmc_machine_signature()
try
    ncores = feature('numcores');
catch
    ncores = 1;
end
sig = sprintf('%s|%d|%s', computer('arch'), ncores, computer);
end

% ------------------------------------------------------------------------
function [n, bytes, secs] = ctmc_profile_point(g, bytesPerNz)
% Factorize the nonsingular block of a g-by-g nearest-neighbour lattice
% generator (n=g*g states). This QBD-like structure is representative of
% multi-station queueing generators and cheap and safe to factorize.
A = ctmc_lattice_generator(g);
n = size(A, 1);
t0 = tic;
[L, U] = lu(A);
secs = toc(t0);
bytes = bytesPerNz * (nnz(L) + nnz(U));
end

% ------------------------------------------------------------------------
function A = ctmc_lattice_generator(g)
n = g * g;
idx = reshape(1:n, g, g);
src = []; dst = [];
% right neighbour
s = idx(:,1:end-1); d = idx(:,2:end);
src = [src; s(:)]; dst = [dst; d(:)];
% left neighbour
s = idx(:,2:end); d = idx(:,1:end-1);
src = [src; s(:)]; dst = [dst; d(:)];
% down neighbour
s = idx(1:end-1,:); d = idx(2:end,:);
src = [src; s(:)]; dst = [dst; d(:)];
% up neighbour
s = idx(2:end,:); d = idx(1:end-1,:);
src = [src; s(:)]; dst = [dst; d(:)];
Q = sparse(src, dst, 1, n, n);
rowsum = full(sum(Q, 2));
Q = Q - spdiags(rowsum, 0, n, n);
A = Q(1:n-1, 1:n-1);   % drop last state -> nonsingular
end

% ------------------------------------------------------------------------
function [alpha, beta] = ctmc_fit_power_law(x1, y1, x2, y2)
if x1 <= 0 || x2 <= 0 || y1 <= 0 || y2 <= 0 || x1 == x2
    error('ctmc_memory_gate:degenerateFit', 'degenerate power-law fit');
end
beta = log(y2 / y1) / log(x2 / x1);
alpha = y1 / (x1 ^ beta);
end
