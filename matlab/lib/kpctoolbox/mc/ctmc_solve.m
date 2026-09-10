%{ @file ctmc_solve.m
 %  @brief Equilibrium distribution of the continuous-time Markov chain
 %
 %  @author LINE Development Team
%}

%{
 % @brief Equilibrium distribution of the continuous-time Markov chain
 %
 % @details
 % Calculates the equilibrium distribution of a continuous-time Markov chain given its infinitesimal generator matrix.
 %
 % @par Syntax:
 % @code
 % p = ctmc_solve(Q)
 % [p, Q, nConnComp, connComp] = ctmc_solve(Q, options)
 % @endcode
 %
 % @par Parameters:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>Q<td>Infinitesimal generator matrix of the continuous-time Markov chain
 % <tr><td>options<td>(Optional) Solver options (config.linsolver, or method as a fallback: 'gmres', 'bicgstab', 'direct', 'gpu' or default; force: boolean, verbose: 2 for debug)
 % </table>
 %
 % @par Returns:
 % <table>
 % <tr><th>Name<th>Description
 % <tr><td>p<td>Equilibrium distribution vector
 % <tr><td>Q<td>Processed generator matrix (e.g., after removing spurious zeros)
 % <tr><td>nConnComp<td>Number of connected components found (if reducible)
 % <tr><td>connComp<td>Vector assigning each state to a connected component
 % </table>
 %
 % @par Examples:
 % @code
 % Q = [-0.5, 0.5; 0.2, -0.2];
 % p = ctmc_solve(Q);
 % @endcode
%}
function [p, Q, nConnComp, connComp]=ctmc_solve(Q,options)

% Order above which the direct sparse factorization is abandoned in favour of the
% Krylov path (GMRES, then BiCGSTAB). The former blocking prompt at this size is
% gone: it warned before a solve that would exhaust memory, and there is now an
% iterative path that does not, with the direct solve retained as the fallback
% when both iterative methods fail.
GMRES_MIN_STATES = 6000;

if size(Q)==1
    p = 1;
    nConnComp = 1;
    connComp = 1:length(Q);
    return
end

Q = ctmc_makeinfgen(Q); % so that spurious diagonal elements are set to 0
n = length(Q);

if issym(Q) && nargin > 1 && isfield(options,'config') && isfield(options.config,'symbolic') ...
        && (strcmpi(options.config.symbolic,'sage') || strncmpi(options.config.symbolic,'http',4))
    % Symbolic solve delegated to the computer algebra backend (SAGE.m). The
    % same request from MATLAB, the JAR and native Python then returns the
    % same normal form, which is what makes symbolic results comparable
    % across the three codebases. The toolbox path below is unchanged and
    % stays the default.
    [p, ~, ~, nConnComp, connComp] = SAGE.solveCTMC(Q, {}, ...
        SAGE.resolve(options.config.symbolic));
    p = reshape(p, 1, []);
    return
end

if issym(Q)
    symvariables = symvar(Q); % find all symbolic variables
    B = double(subs(Q+Q',symvariables,ones(size(symvariables)))); % replace all symbolic variables with 1.0
else
    B = abs(Q+Q')>0;
end
[nConnComp, connComp] = weaklyconncomp(B);
if nConnComp > 1
    % reducible generator - solve each component recursively
    line_warning(mfilename,'Reducible generator. No initial vector available, decomposing and solving each component recursively.\n');
    if issym(Q)
        p = sym(zeros(1,n));
    else
        p = zeros(1,n);
    end

    for c=1:nConnComp
        Qc = Q(connComp==c,connComp==c);
        Qc = ctmc_makeinfgen(Qc);
        p(connComp==c) = ctmc_solve(Qc);
    end
    p = p /sum(p);
    return
end

if all(Q==0)
    % No transitions at all: every distribution satisfies p*Q=0, so the
    % stationary distribution is not unique and uniform is as good as any.
    p = ones(1,n)/n;
    return
end
p = zeros(1,n);
b = zeros(n,1);

nnzel = 1:n;
Qnnz = Q; bnnz = b;
Qnnz_1 = Qnnz; bnnz_1 = bnnz;

isReducible = false;
goon = true;
% AN ISOLATED STATE IS DROPPED; AN ABSORBING ONE IS NOT. The column mass
% sum(abs(Q(:,j))) counts the inflow of j plus its own outflow through the
% diagonal, so it vanishes exactly for a state with neither -- an isolated
% state, which carries no stationary mass and only makes the system singular.
% The test used to ALSO require a nonzero ROW mass, which vanishes exactly for
% an ABSORBING state: the one state the stationary mass ends up in. Dropping it
% left the states that fed it with nothing to flow into, CTMC_MAKEINFGEN then
% re-zeroed their diagonals, and the elimination cascaded until nothing was
% left and this function raised "no recurrent state" on a chain whose
% stationary distribution is perfectly unique (Q = [0 0; 1 -1] has pi = [1 0]).
%
% THAT COST A HOST-DEPENDENT ANSWER, not just a refusal. The row mass is a
% COMPUTED SUM tested against an exact zero, so a generator assembled slightly
% differently on two CPUs -- the same MAP built through a different BLAS kernel
% -- lands on 0 for one and 1e-17 for the other, and the two then take opposite
% branches. `dec.source.mmap` on the self-looping sanity models is where this
% surfaced: MAM's traffic merge asks MMAP_LAMBDA for the rate of a link whose
% phase process settles in one phase, and the refusal reached the bisection in
% SOLVER_MAM_BASIC_MMAP_CLOSED as "that arrival rate overloads the network",
% which drove it to a different lambda on picard04 (Tput 0.2296) than on
% picard09 (0.6732) for the same model and code. See _kb/06-solver-catalog.md.
%
% Native python already draws the line here -- `col_sums < 1e-12` in
% `api/mc/ctmc.py`, with the comment "Not the same as absorbing, which has
% incoming transitions and a zero off-diagonal ROW" -- so this is MATLAB
% catching up rather than a new convention.
isolatedTol = 1e-12;
while goon
    colmass = sum(abs(Qnnz),1);
    if issym(Qnnz)
        nnzel = find(colmass ~= 0);
    else
        nnzel = find(colmass > isolatedTol);
    end
    if length(nnzel) < n && ~isReducible
        isReducible = true;
        if (nargin > 1 && options.verbose == 2) % debug
            line_warning(mfilename,'The infinitesimal generator is reducible.\n');
        end
    end
    Qnnz = Qnnz(nnzel, nnzel);
    bnnz = bnnz(nnzel);
    Qnnz = ctmc_makeinfgen(Qnnz);
    if all(size(Qnnz_1(:)) == size(Qnnz(:))) && all(size(bnnz_1(:)) == size(bnnz(:)))
        goon = false;
    else
        Qnnz_1 = Qnnz; bnnz_1 = bnnz; nnzel = 1:length(Qnnz);
    end
end

if isempty(Qnnz)
    % Every state was ISOLATED -- no inflow and no outflow anywhere -- so the
    % elimination above emptied the generator. An all-zero Q is answered
    % uniformly further up, so reaching here means the states carried mass in
    % Q but none of it connected. Returning a uniform vector would NOT satisfy
    % p*Q=0 (it is only a shape of the right size) and a caller cannot tell it
    % apart from a real answer: a generator missing all its arrivals reads back
    % as a plausible mean of cutoff/2. Fail instead. A chain with SEVERAL
    % recurrent classes has no unique stationary distribution without an initial
    % vector either, and belongs in ctmc_solve_reducible(Q, pi0); a chain with
    % ONE absorbing state is no longer refused here, its distribution being the
    % point mass the elimination used to throw away.
    line_error(mfilename, sprintf(['The infinitesimal generator has no connected state: every state was eliminated as isolated.\n' ...
        'This generator admits no unique stationary distribution. It usually means the generator is malformed -- ' ...
        'e.g. states that carry a diagonal but no transition between them, as happens when a class of ' ...
        'transitions was dropped while building it. Use ctmc_solve_reducible(Q, pi0) for a chain with several recurrent classes.']));
end
Qnnz_1 = Qnnz;
Qnnz(:,end) = 1;
bnnz_1 = Qnnz;
bnnz(end) = 1;

if ~isdeployed
    if issym(Q)
        p = sym(p);
    end
end

warning('off','MATLAB:singularMatrix');

% Iterative path. The direct solve stays the default and remains the last
% fallback: a Krylov method is used only above GMRES_MIN_STATES, or when
% explicitly requested, and only when it reports convergence. A symbolic
% generator always takes the direct path, there being no iterative method over a
% symbolic field.
%
% Two Krylov methods are tried in sequence before the direct solve. GMRES(m) is
% first because its residual is monotone and it is the more robust of the two.
% BiCGSTAB follows because the way GMRES(m) fails on a generator is stagnation,
% the useful subspace being wider than the restart window, and a short-recurrence
% method has no restart to stagnate on. The direct solve is cubic at this size,
% so a second iterative attempt is cheap against what it may avoid.
% WHICH LINEAR SOLVE, not which method. 'gmres', 'bicgstab', 'direct' and
% 'gpu' pick a backend for the SAME generator, so they are read from
% OPTIONS.CONFIG.LINSOLVER, which is where a LINE solver can set them:
% SolverCTMC's method namespace is about state-space construction (default,
% exact, gpu, mdd, cftp) and runAnalyzerChecks refuses a name outside it, so a
% backend named through options.method could never reach here from a solver.
% OPTIONS.METHOD is still honoured as the fallback: this is a kpctoolbox
% library function and its other callers (mdd_closedqn's ctmcmethod, direct
% callers passing struct('method',...)) name the backend that way.
method = 'default';
if nargin > 1 && isfield(options,'config') && isstruct(options.config) ...
        && isfield(options.config,'linsolver') && ~isempty(options.config.linsolver)
    method = lower(char(options.config.linsolver));
elseif nargin > 1 && isfield(options,'method') && ~isempty(options.method)
    method = lower(options.method);
end
isKrylovName = strcmp(method,'gmres') || strcmp(method,'bicgstab');
useKrylov = ~issym(Q) && (isKrylovName || ...
    (~strcmp(method,'direct') && length(Qnnz) > GMRES_MIN_STATES));
if useKrylov
    restart = [];
    if nargin > 1 && isfield(options,'config') && isfield(options.config,'gmres_restart')
        restart = options.config.gmres_restart;
    end
    maxit = [];
    if nargin > 1 && isfield(options,'iter_max') && ~isempty(options.iter_max)
        if isempty(restart)
            maxit = min(ceil(length(Qnnz)/min(length(Qnnz),50)), options.iter_max);
        else
            maxit = min(ceil(length(Qnnz)/restart), options.iter_max);
        end
    end
    verbose2 = nargin > 1 && isfield(options,'verbose') && options.verbose == 2;
    if ~strcmp(method,'bicgstab')
        [xg,gflag] = ctmc_gmres(Qnnz', bnnz, [], restart, maxit, []);
        if gflag == 0
            p(nnzel) = xg;
            warning('on','MATLAB:singularMatrix');
            return
        end
        if verbose2
            line_warning(mfilename,'GMRES did not converge (flag %d), trying BiCGSTAB.\n', gflag);
        end
    end
    bmaxit = [];
    if nargin > 1 && isfield(options,'iter_max') && ~isempty(options.iter_max)
        bmaxit = options.iter_max;
    end
    [xb,bflag] = ctmc_bicgstab(Qnnz', bnnz, [], bmaxit, []);
    if bflag == 0
        p(nnzel) = xb;
        warning('on','MATLAB:singularMatrix');
        return
    end
    if verbose2
        line_warning(mfilename,'BiCGSTAB did not converge (flag %d), falling back to the direct solve.\n', bflag);
    end
end

if nargin == 1
    p(nnzel)=Qnnz'\ bnnz;
    if any(isnan(p))
        % verify if this has become reducible        
        if issym(Qnnz)
            symvariables = symvar(Qnnz); % find all symbolic variables
            B = double(subs(Qnnz+Qnnz',symvariables,ones(size(symvariables)))); % replace all symbolic variables with 1.0
        else
            B = abs(Qnnz+Qnnz')>0;
        end
        [nConnComp, connComp] = weaklyconncomp(B);
        if nConnComp > 1
            % reducible generator - solve each component recursively
            if issym(Qnnz)
                p(nnzel) = sym(zeros(1,n));
            else
                p(nnzel) = zeros(1,n);
            end

            for c=1:nConnComp
                Qc = Q(connComp==c,connComp==c);
                Qc = ctmc_makeinfgen(Qc);
                p(intersect(find(connComp==c),nnzel)) = ctmc_solve(Qc);
            end
            p = p /sum(p);            
            return
        end
        % ONE COMPONENT AND STILL SINGULAR MEANS SEVERAL RECURRENT CLASSES.
        % pi*Q = 0 then has a solution SPACE rather than a solution, and which
        % member the factorization lands on says nothing about the chain -- it
        % depends on the initial distribution, which this signature does not
        % carry. Returning the NaNs would push the ambiguity into the caller's
        % arithmetic silently, so refuse here as this function did before
        % absorbing states were kept.
        line_error(mfilename, sprintf(['The infinitesimal generator admits no unique stationary distribution: '...
            'the balance equations are singular over a single connected component, which is what several '...
            'recurrent classes look like. Use ctmc_solve_reducible(Q, pi0), which resolves the ambiguity '...
            'with an initial vector.']));
    end
else
    % same backend selector resolved above (config.linsolver, else method)
    switch method
        case 'gpu'
            try
                gQnnz = gpuArray(Qnnz');
                gbnnz = gpuArray(bnnz);
                pGPU = gQnnz \ gbnnz;
                gathered_pGPU = gather(pGPU);
                p(nnzel) = gathered_pGPU; % transfer from GPU to local env
            catch
                warning('ctmc_solve: GPU either not available or execution failed. Switching to default method.');
                p(nnzel) = Qnnz'\ bnnz;
            end
        otherwise
            p(nnzel)=Qnnz'\ bnnz;
    end
end

if issym(Q)
    Q=simplify(Q);
end
warning('on','MATLAB:singularMatrix');
end
