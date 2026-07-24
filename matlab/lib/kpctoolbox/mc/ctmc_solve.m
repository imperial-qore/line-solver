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
 % <tr><td>options<td>(Optional) Solver options (method: 'gpu' or default, force: boolean, verbose: 2 for debug)
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

% Order above which the direct sparse factorization is abandoned in favour of
% GMRES. The former blocking prompt at this size is gone: it warned before a
% solve that would exhaust memory, and there is now an iterative path that does
% not, with the direct solve retained as the fallback when GMRES fails.
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
while goon
    nnzel = find(sum(abs(Qnnz),1)~=0 & sum(abs(Qnnz),2)'~=0);
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
    % The elimination above drops every state whose row is all-zero, which is
    % precisely an ABSORBING state; ctmc_makeinfgen then re-zeroes the diagonal
    % of the survivors that only fed it, so the elimination cascades until
    % nothing is left. Returning a uniform vector here does NOT satisfy p*Q=0
    % (it is not a stationary distribution, just a shape of the right size), and
    % a caller cannot tell it apart from a real answer: a generator missing all
    % its arrivals reads back as a plausible mean of cutoff/2. Fail instead.
    % A genuinely absorbing chain has no unique stationary distribution without
    % an initial vector, so it belongs in ctmc_solve_reducible(Q, pi0).
    line_error(mfilename, sprintf(['The infinitesimal generator has no recurrent state: every state was eliminated as absorbing.\n' ...
        'This generator admits no unique stationary distribution. It usually means the generator is malformed -- ' ...
        'e.g. a state with no outgoing transitions that absorbs the whole chain, as happens when a class of ' ...
        'transitions was dropped while building it. Use ctmc_solve_reducible(Q, pi0) for a genuinely absorbing chain.']));
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

% Iterative path. The direct solve stays the default and remains the fallback:
% GMRES is used only above GMRES_MIN_STATES, or when explicitly requested, and
% only when it reports convergence. A symbolic generator always takes the direct
% path, there being no iterative method over a symbolic field.
method = 'default';
if nargin > 1 && isfield(options,'method') && ~isempty(options.method)
    method = lower(options.method);
end
useGmres = ~issym(Q) && (strcmp(method,'gmres') || ...
    (~strcmp(method,'direct') && length(Qnnz) > GMRES_MIN_STATES));
if useGmres
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
    [xg,gflag] = ctmc_gmres(Qnnz', bnnz, [], restart, maxit, []);
    if gflag == 0
        p(nnzel) = xg;
        warning('on','MATLAB:singularMatrix');
        return
    end
    if nargin > 1 && isfield(options,'verbose') && options.verbose == 2
        line_warning(mfilename,'GMRES did not converge (flag %d), falling back to the direct solve.\n', gflag);
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
    end
else
    if ~isfield(options, 'method')
        options.method = 'default';
    end
    switch options.method
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
