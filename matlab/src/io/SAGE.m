classdef SAGE
    % SAGE  MATLAB-to-SageMath bridge for symbolic analysis.
    %
    % Static helpers that talk to the line-sage-rest service, the same JSON
    % protocol the JAR (jline.api.sym) and native Python
    % (line_solver.api.sym) use. The service is SageMath in a container; see
    % io/sage/server.py for the endpoints.
    %
    % It exists for two reasons. The Symbolic Math Toolbox is licence-gated,
    % so a session without it cannot run any symbolic analysis at all; and
    % where both are present, routing through one engine gives all three
    % codebases the same normal form, which is what makes a symbolic result
    % comparable across them.
    %
    % Backend resolution, in order:
    %   1. an explicit URL in options.config.symbolic;
    %   2. the LINE_SAGE_URL environment variable;
    %   3. a line-sage-rest service already listening on a conventional port;
    %   4. a container started here from a locally present image;
    %   5. nothing, in which case the caller falls back to the toolbox.
    %
    % Step 3 checks identity through /api/v1/info rather than trusting the
    % port: every imperialqore line-*-rest service listens on 8080 by
    % convention, so a health probe alone would accept the LQNS service.
    %
    % Expressions cross as plain char infix strings, e.g. '2*x1 - 3*x2'.
    % Numeric coefficients are read server side as exact rationals, so
    % '0.1' is 1/10 and not the binary double nearest to it.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        DOCKER_IMAGES = {'imperialqore/line-sage-rest:latest', ...
            'imperialqore/line-sage-rest'};
        PROBE_PORTS = [8085, 8080];
        STARTUP_TIMEOUT = 120;   % seconds to wait for a container to answer
        DEFAULT_TIMEOUT = 300;   % seconds per request
        % The POST routes this client calls; a service that does not serve
        % them all is older than the client and is refused by name. See
        % servesRequiredRoutes.
        REQUIRED_ROUTES = {'/api/v1/ctmc/solve', '/api/v1/ctmc/sensitivity', ...
            '/api/v1/ctmc/measures', '/api/v1/ctmc/passage', '/api/v1/simplify', ...
            '/api/v1/diff', '/api/v1/eval', '/api/v1/fluid/odes'};
        % A service that answers /health can still be unable to COMPUTE. The
        % line-sage-rest image ships a FLINT built for CPUs that have BMI2 and
        % ADX; on an older host the first multi-limb exact operation raises
        % SIGILL, the worker dies mid-request and the call returns no bytes at
        % all. Health is pure Python and keeps answering, so it cannot see
        % this. The canary is the weighted-average softmin form -- what the
        % fluid export actually sends -- and is the smallest expression
        % observed to trigger it; see _kb/11-conventions-and-gotchas.md.
        CANARY_EXPR = '(x*exp(-x) + exp(-1))/(exp(-x) + exp(-1))';
        CANARY_ARG = '0.68999999999999995';
        CANARY_VALUE = 0.82116556904906557;
    end

    methods (Static)

        function url = resolve(requested)
            % RESOLVE  Base URL of a usable service, or '' if there is none.
            %
            % REQUESTED is '' or 'auto' to search, a URL to use a specific
            % service, 'none' to disable the backend, or an image name.
            persistent cachedUrl
            if nargin < 1 || isempty(requested)
                requested = 'auto';
            end
            requested = strtrim(requested);
            if any(strcmpi(requested, {'none', 'off'}))
                url = '';
                return
            end
            if strncmpi(requested, 'http://', 7) || strncmpi(requested, 'https://', 8)
                url = requested;
                if ~SAGE.isReachable(url) || ~SAGE.servesRequiredRoutes(url) || ~SAGE.isUsable(url)
                    url = '';
                end
                return
            end

            env = getenv('LINE_SAGE_URL');
            if ~isempty(env) && SAGE.isReachable(env) && SAGE.servesRequiredRoutes(env) ...
                    && SAGE.isUsable(env)
                url = env;
                return
            end

            if ~isempty(cachedUrl) && SAGE.isReachable(cachedUrl) ...
                    && SAGE.servesRequiredRoutes(cachedUrl) && SAGE.isUsable(cachedUrl)
                url = cachedUrl;
                return
            end

            for h = SAGE.probeHosts()
                for p = SAGE.PROBE_PORTS
                    candidate = sprintf('http://%s:%d', h{1}, p);
                    if SAGE.isSageService(candidate) && SAGE.servesRequiredRoutes(candidate) ...
                            && SAGE.isUsable(candidate)
                        url = candidate;
                        cachedUrl = url;
                        return
                    end
                end
            end

            if any(strcmpi(requested, {'auto', 'true', 'sage'}))
                % 'sage' names the ENGINE, not a Docker image: taking it as an
                % image name sends docker looking for a repository called
                % "sage" and the start fails with a pull-access error that
                % reads like a login problem.
                image = SAGE.getDockerImage();
            else
                image = requested;
            end
            % Auto-pull only on an explicit opt-in: the 'sage' keyword or a named
            % image. Bare 'auto'/'true'/'' keep the native backend unless the
            % image is already local, so leaving symbolic on auto never pulls.
            if isempty(image) && strcmpi(requested, 'sage')
                image = SAGE.pullDockerImage(SAGE.DOCKER_IMAGES{1});
            elseif ~isempty(image) && ~any(strcmpi(requested, {'auto', 'true', 'sage'})) ...
                    && ~SAGE.hasLocalImage(image)
                image = SAGE.pullDockerImage(image);
            end
            if isempty(image)
                url = '';
                return
            end
            url = SAGE.startContainer(image);
            cachedUrl = url;
        end

        function hosts = probeHosts()
            % PROBEHOSTS  Host names to search for a service, in order.
            %
            % INSIDE THE CONTAINER, `localhost` IS THE CONTAINER. matlab_docker.sh
            % gives MATLAB its own network namespace on purpose (the MathWorks
            % Service Host claims an abstract AF_UNIX socket, which is scoped to
            % that namespace), so a line-sage-rest published on the host's port
            % answers on host.docker.internal and on nothing else this process
            % can name. Without this the whole symbolic ladder skipped under
            % test phase 6 on hosts that had the image and could have run it,
            % and startContainer could not help: there is no docker in here.
            % The alias comes from the --add-host in matlab_docker.sh; an older
            % container without it simply fails the probe, as before.
            hosts = {'localhost'};
            if ~isempty(getenv('LINE_MATLAB_IN_CONTAINER'))
                hosts{end+1} = 'host.docker.internal';
            end
        end

        function bool = isAvailable(requested)
            % ISAVAILABLE  True if a symbolic backend can be resolved.
            if nargin < 1
                requested = 'auto';
            end
            bool = ~isempty(SAGE.resolve(requested));
        end

        function img = getDockerImage()
            % GETDOCKERIMAGE  First locally present image tag, or '' if none.
            img = '';
            if ispc
                return
            end
            if unix('docker info >/dev/null 2>&1') ~= 0
                return
            end
            for i = 1:numel(SAGE.DOCKER_IMAGES)
                [st, out] = unix(['docker images -q ', SAGE.DOCKER_IMAGES{i}, ' 2>/dev/null']);
                if st == 0 && ~isempty(strtrim(out))
                    img = SAGE.DOCKER_IMAGES{i};
                    return
                end
            end
        end

        function tf = hasLocalImage(image)
            % HASLOCALIMAGE  True if IMAGE is present in the local Docker store.
            tf = false;
            if ispc || isempty(image)
                return
            end
            [st, out] = unix(['docker images -q ', image, ' 2>/dev/null']);
            tf = (st == 0) && ~isempty(strtrim(out));
        end

        function img = pullDockerImage(target)
            % PULLDOCKERIMAGE  Pull TARGET if Docker is usable and there is room.
            %
            % Returns the tag on success, else ''. Storage-guarded via
            % lineDockerHasStorageFor (the same guard as the LQNS/JMT
            % wrappers), so an opt-in symbolic request never silently fills the
            % Docker disk; on refusal the caller keeps its native algebra.
            img = '';
            if ispc
                return
            end
            if unix('docker info >/dev/null 2>&1') ~= 0
                return
            end
            if ~lineDockerHasStorageFor(target)
                line_warning(mfilename, ['Skipping docker pull of %s: insufficient free ' ...
                    'space at the Docker storage location; keeping the native symbolic ' ...
                    'backend.\n'], target);
                return
            end
            line_printf('[LINE] Pulling Docker image %s (this may take a while)...\n', target);
            if unix(['docker pull ', target]) == 0 && SAGE.hasLocalImage(target)
                img = target;
            end
        end

        function out = startedContainer(name)
            % STARTEDCONTAINER  Get/set the container this session started.
            %
            % Call with a name to set, '' to clear, or no argument to get the
            % tracked name (parity with the Java/Python startedContainer state).
            persistent tracked
            if nargin >= 1
                tracked = name;
            end
            if isempty(tracked)
                out = '';
            else
                out = tracked;
            end
        end

        function url = startContainer(image)
            % STARTCONTAINER  Run the service and wait for it to answer.
            %
            % The container is bound to a free host port, so several MATLAB
            % sessions, or a session next to a hand-started service, do not
            % collide. It stays warm for the session (repeated symbolic calls
            % are then cheap) and is best-effort stopped when MATLAB exits or
            % SAGE is cleared, via the onCleanup guard below. Stop it sooner with
            % SAGE.stopContainer.
            % Held only for its destructor: cleared at MATLAB exit / clear SAGE,
            % which fires the onCleanup below and stops the container.
            persistent guard
            url = '';
            port = SAGE.freePort();
            name = sprintf('line-sage-rest-%d', port);
            cmd = sprintf('docker run -d --rm --name %s -p %d:8080 %s', name, port, image);
            [st, out] = unix([cmd, ' 2>&1']);
            if st ~= 0
                line_warning(mfilename, 'Could not start %s: %s\n', image, strtrim(out));
                return
            end
            % Track only this container so stopContainer targets it precisely
            % (not every line-sage-rest-* on the host), matching Java/Python.
            SAGE.startedContainer(name);
            guard = onCleanup(@() SAGE.stopContainerByName(name)); %#ok<NASGU>
            candidate = sprintf('http://localhost:%d', port);
            deadline = tic;
            while toc(deadline) < SAGE.STARTUP_TIMEOUT
                if SAGE.isReachable(candidate)
                    if SAGE.servesRequiredRoutes(candidate) && SAGE.isUsable(candidate)
                        url = candidate;
                        return
                    end
                    % Booted, but its arithmetic dies on this CPU. Keeping it
                    % running would only cost memory, and returning it would
                    % hand the caller a backend that kills every request.
                    SAGE.stopContainerByName(name);
                    SAGE.startedContainer('');
                    return
                end
                pause(0.5);
            end
            SAGE.stopContainerByName(name);
            SAGE.startedContainer('');
            line_warning(mfilename, ...
                'Container %s did not answer within %d s.\n', name, SAGE.STARTUP_TIMEOUT);
        end

        function stopContainer()
            % STOPCONTAINER  Stop the container this session started, if any.
            %
            % Parity with the Java (shutdown hook) and Python (atexit) clients:
            % it stops only the container started here, not every line-sage-rest-*
            % on the host.
            name = SAGE.startedContainer();
            if ~isempty(name)
                SAGE.stopContainerByName(name);
                SAGE.startedContainer('');
            end
        end

        function stopContainerByName(name)
            % STOPCONTAINERBYNAME  Stop a single named container, best effort.
            if ispc || isempty(name)
                return
            end
            unix(sprintf('docker stop -t 1 %s >/dev/null 2>&1', name));
        end

        function bool = isReachable(url)
            % ISREACHABLE  True if the service answers a health probe.
            bool = false;
            try
                opts = weboptions('Timeout', 5, 'ContentType', 'json');
                r = webread([SAGE.trimUrl(url), '/api/v1/health'], opts);
                bool = isfield(r, 'status') && strcmp(r.status, 'ok');
            catch
            end
        end

        function bool = isUsable(url)
            % ISUSABLE  True if the service answers AND can evaluate.
            %
            % Verdicts are cached per URL: this costs one small request the
            % first time a service is considered, and nothing after.
            persistent verdicts
            if isempty(verdicts)
                verdicts = containers.Map('KeyType', 'char', 'ValueType', 'logical');
            end
            key = SAGE.trimUrl(url);
            if isKey(verdicts, key)
                bool = verdicts(key);
                return
            end
            bool = false;
            try
                opts = weboptions('MediaType', 'application/json', ...
                    'ContentType', 'json', 'Timeout', 60);
                payload = struct('exprs', {{SAGE.CANARY_EXPR}}, ...
                    'values', struct('x', SAGE.CANARY_ARG), 'timeout_s', 30);
                r = webwrite([key, '/api/v1/eval'], payload, opts);
                if isfield(r, 'status') && strcmp(r.status, 'ok') && isfield(r, 'values')
                    v = r.values;
                    if iscell(v)
                        v = v{1};
                    else
                        v = v(1);
                    end
                    bool = isnumeric(v) && isfinite(v) ...
                        && abs(v - SAGE.CANARY_VALUE) < 1e-9;
                end
            catch
                % A dead worker closes the connection without a reply, which
                % surfaces here as a webservices error rather than a service
                % one. Either way the backend cannot serve us.
            end
            if ~bool
                line_warning(mfilename, ...
                    ['Ignoring symbolic backend at %s: it did not return the ', ...
                     'usability canary. On a CPU without BMI2/ADX the image''s ', ...
                     'FLINT raises SIGILL mid-request.\n'], key);
            end
            verdicts(key) = bool;
        end

        function bool = isSageService(url)
            % ISSAGESERVICE  True if the service is line-sage-rest and not
            % another line-*-rest service sharing the conventional port.
            bool = false;
            try
                opts = weboptions('Timeout', 5, 'ContentType', 'json');
                r = webread([SAGE.trimUrl(url), '/api/v1/info'], opts);
                bool = isfield(r, 'sage_version');
            catch
            end
        end

        function bool = servesRequiredRoutes(url)
            % SERVESREQUIREDROUTES  True if the service is new enough to use.
            %
            % A SERVICE CAN BE HEALTHY AND STILL TOO OLD, and that case used to
            % reach the caller as a request failure rather than as a resolution
            % failure. The roster ran a line-sage-rest image predating the
            % ctmc/measures and ctmc/passage routes on 2026-09-11: it answered
            % /api/v1/health, passed the isUsable arithmetic canary, and then
            % failed every call to those routes with "unknown endpoint
            % /api/v1/ctmc/measures" -- an error that reads as a defect in the
            % solver rather than as an out-of-date container. /api/v1/info
            % lists what the running server actually routes, so asking it turns
            % that into a named refusal here, where the caller falls back to
            % its own algebra or skips and is told why.
            bool = false;
            try
                opts = weboptions('Timeout', 5, 'ContentType', 'json');
                r = webread([SAGE.trimUrl(url), '/api/v1/info'], opts);
            catch
                return
            end
            if ~isfield(r, 'endpoints')
                % A service that does not list its routes cannot be
                % interrogated; leave it to the request itself to fail.
                bool = true;
                return
            end
            served = r.endpoints;
            if ~iscell(served)
                served = cellstr(served);
            end
            missing = SAGE.REQUIRED_ROUTES(~ismember(SAGE.REQUIRED_ROUTES, served(:)'));
            if isempty(missing)
                bool = true;
                return
            end
            line_warning(mfilename, sprintf(['Ignoring the line-sage-rest service at %s: it is ', ...
                'older than this client and does not serve %s. Refresh it with ', ...
                'docker pull %s.'], SAGE.trimUrl(url), strjoin(missing, ', '), ...
                SAGE.DOCKER_IMAGES{1}));
        end

        % ---------------------------------------------------------------
        % Operations
        % ---------------------------------------------------------------

        function [pi, num, den, nConnComp, connComp] = solveCTMC(Q, symbols, url, timeout)
            % SOLVECTMC  Symbolic stationary distribution, pi*Q = 0, sum(pi) = 1.
            %
            % Q is a cell array of expression strings or a sym matrix, SYMBOLS
            % a cellstr of the symbols occurring in it. PI comes back as a
            % cellstr, or as a sym vector when the Symbolic Toolbox is present.
            if nargin < 3 || isempty(url)
                url = SAGE.require();
            end
            if nargin < 4
                timeout = SAGE.DEFAULT_TIMEOUT;
            end
            [Qcell, symbols] = SAGE.toExpressionMatrix(Q, symbols);
            payload = struct('Q', {Qcell}, 'symbols', {symbols(:)'}, ...
                'normalize', true, 'timeout_s', timeout);
            r = SAGE.post(url, '/api/v1/ctmc/solve', payload, timeout);
            pi = SAGE.toSymOrChar(r.pi);
            num = SAGE.toSymOrChar(r.num);
            den = SAGE.toSymOrChar({r.den});
            if iscell(den)
                den = den{1};
            end
            nConnComp = double(r.nConnComp);
            connComp = double(r.connComp(:))';
        end

        function [dpi, S, SS, pi] = sensitivity(Q, symbols, theta, reward, url, timeout)
            % SENSITIVITY  Exact d(pi)/d(theta) and, given a reward, d(E[r])/d(theta).
            %
            % Exact where @SolverCTMC/getSensitivity uses a central difference
            % accurate to O(h^2). As there, dr/dtheta is taken to be zero.
            if nargin < 5 || isempty(url)
                url = SAGE.require();
            end
            if nargin < 6
                timeout = SAGE.DEFAULT_TIMEOUT;
            end
            [Qcell, symbols] = SAGE.toExpressionMatrix(Q, symbols);
            payload = struct('Q', {Qcell}, 'symbols', {symbols(:)'}, ...
                'theta', theta, 'timeout_s', timeout);
            if nargin >= 4 && ~isempty(reward)
                payload.reward = SAGE.toExpressionList(reward);
            end
            r = SAGE.post(url, '/api/v1/ctmc/sensitivity', payload, timeout);
            dpi = SAGE.toSymOrChar(r.dpi);
            pi = SAGE.toSymOrChar(r.pi);
            S = '';
            SS = '';
            if isfield(r, 'S')
                S = SAGE.scalarSymOrChar(r.S);
            end
            if isfield(r, 'SS')
                SS = SAGE.scalarSymOrChar(r.SS);
            end
        end

        function [meas, ratios, pi, num, den] = ctmcMeasures(Q, symbols, weights, names, ratios, url, timeout)
            % CTMCMEASURES  Stationary distribution and weighted sums of it.
            %
            % [MEAS, RATIOS, PI, NUM, DEN] = CTMCMEASURES(Q, SYMBOLS, WEIGHTS,
            %   NAMES, RATIOS, URL, TIMEOUT)
            %
            % Every mean SolverCTMC reports is a NUMERIC linear functional of
            % pi, so the caller builds one weight vector per measure with no
            % algebra at all and this single round trip returns each w*pi as a
            % rational function: a throughput is pi*depRates, a queue length is
            % pi*stateSpaceAggr, a marginal probability is pi*indicator and a
            % mean reward is pi*r.
            %
            % WEIGHTS is (nmeasures x nstates), numeric or a cellstr of
            % expressions: the entries are parsed in the same field as Q, so an
            % expression is admissible and not only a number. RATIOS is a struct
            % array with fields name, num and den, the last two being 0-based
            % indices into WEIGHTS. A ratio whose denominator measure is
            % identically zero comes back with expr 'undefined' and a reason,
            % never as 0.
            %
            % MEAS and RATIOS are ALWAYS struct arrays with the fields name,
            % expr, num, den and reason, '' where the service sent nothing. The
            % normalization is not cosmetic: a defined measure and an undefined
            % one carry different JSON fields, and jsondecode answers a cell
            % array rather than a struct array whenever they are mixed.
            if nargin < 6 || isempty(url)
                url = SAGE.require();
            end
            if nargin < 7
                timeout = SAGE.DEFAULT_TIMEOUT;
            end
            [Qcell, symbols] = SAGE.toExpressionMatrix(Q, symbols);
            payload = struct('Q', {Qcell}, 'symbols', {symbols(:)'}, ...
                'weights', {SAGE.toExpressionRows(weights)}, ...
                'normalize', true, 'timeout_s', timeout);
            if nargin >= 4 && ~isempty(names)
                payload.names = SAGE.asCellstr(names);
            end
            if nargin >= 5 && ~isempty(ratios)
                payload.ratios = ratios;
            end
            r = SAGE.post(url, '/api/v1/ctmc/measures', payload, timeout);
            meas = SAGE.toMeasureStruct(r.measures);
            ratios = SAGE.toMeasureStruct([]);
            if isfield(r, 'ratios')
                ratios = SAGE.toMeasureStruct(r.ratios);
            end
            pi = SAGE.toSymOrChar(r.pi);
            num = SAGE.toSymOrChar(r.num);
            den = SAGE.toSymOrChar({r.den});
            if iscell(den)
                den = den{1};
            end
        end

        function out = ctmcPassage(S, s0, alpha, atom, symbols, svar, want, nmax, url, timeout)
            % CTMCPASSAGE  First passage time into a target set.
            %
            % OUT = CTMCPASSAGE(S, S0, ALPHA, ATOM, SYMBOLS, SVAR, WANT, NMAX,
            %   URL, TIMEOUT)
            %
            % With A the complement of the target, S = Q(A,A) the sub-generator,
            % S0 = -S*1 the exit vector and ALPHA the initial law restricted to
            % A, following Harrison and Knottenbelt (2002),
            %
            %   L(s) = alpha*(s*I - S)^-1*s0 + atom          Eqs. 1-2
            %   (-S)*M(n) = n*M(n-1),  M(0) = 1              Eq. 3
            %
            % The MOMENTS come from the recursion, not from differentiating L:
            % it carries no transform symbol, stays in the same fraction field
            % and is NMAX right solves against one matrix. WANT is any of
            % 'lst', 'lstall', 'moments', 'momall'. SVAR is adjoined to the
            % field only when a transform is asked for and must differ from
            % every rate symbol.
            if nargin < 4 || isempty(atom)
                atom = '0';
            end
            if nargin < 5
                symbols = {};
            end
            if nargin < 6 || isempty(svar)
                svar = 's';
            end
            if nargin < 7 || isempty(want)
                want = {'moments'};
            end
            if nargin < 8 || isempty(nmax)
                nmax = 1;
            end
            if nargin < 9 || isempty(url)
                url = SAGE.require();
            end
            if nargin < 10
                timeout = SAGE.DEFAULT_TIMEOUT;
            end
            [Scell, symbols] = SAGE.toExpressionMatrix(S, symbols);
            payload = struct('S', {Scell}, ...
                's0', {SAGE.toExpressionList(s0)}, ...
                'alpha', {SAGE.toExpressionList(alpha)}, ...
                'atom', SAGE.exprToChar(atom), ...
                'symbols', {symbols(:)'}, 'svar', svar, ...
                'want', {SAGE.asCellstr(want)}, 'nmax', double(nmax), ...
                'timeout_s', timeout);
            out = SAGE.post(url, '/api/v1/ctmc/passage', payload, timeout);
        end

        function out = simplify(exprs, form, url, timeout)
            % SIMPLIFY  Rewrite expressions: simplify, factor, together,
            % cancel, expand or latex.
            if nargin < 2 || isempty(form)
                form = 'cancel';
            end
            if nargin < 3 || isempty(url)
                url = SAGE.require();
            end
            if nargin < 4
                timeout = SAGE.DEFAULT_TIMEOUT;
            end
            payload = struct('exprs', {SAGE.toExpressionList(exprs)}, ...
                'form', form, 'timeout_s', timeout);
            r = SAGE.post(url, '/api/v1/simplify', payload, timeout);
            out = SAGE.asCellstr(r.results);
        end

        function out = diff(exprs, variable, order, url, timeout)
            % DIFF  Differentiate expressions with respect to VARIABLE.
            if nargin < 3 || isempty(order)
                order = 1;
            end
            if nargin < 4 || isempty(url)
                url = SAGE.require();
            end
            if nargin < 5
                timeout = SAGE.DEFAULT_TIMEOUT;
            end
            payload = struct('exprs', {SAGE.toExpressionList(exprs)}, ...
                'var', variable, 'order', order, 'timeout_s', timeout);
            r = SAGE.post(url, '/api/v1/diff', payload, timeout);
            out = SAGE.asCellstr(r.results);
        end

        function [values, exact] = eval(exprs, assignment, url, timeout)
            % EVAL  Substitute values for symbols and evaluate.
            %
            % ASSIGNMENT is a struct whose fields are symbol names. This is how
            % symbolic results are compared across codebases: symbol numbering
            % follows event enumeration order, so comparing expression text is
            % unsound while comparing substituted values is not.
            if nargin < 3 || isempty(url)
                url = SAGE.require();
            end
            if nargin < 4
                timeout = SAGE.DEFAULT_TIMEOUT;
            end
            names = fieldnames(assignment);
            values_in = struct();
            for i = 1:numel(names)
                % Sent as text so the server reads the decimal exactly.
                values_in.(names{i}) = sprintf('%.17g', assignment.(names{i}));
            end
            payload = struct('exprs', {SAGE.toExpressionList(exprs)}, ...
                'values', values_in, 'timeout_s', timeout);
            r = SAGE.post(url, '/api/v1/eval', payload, timeout);
            raw = r.values;
            if iscell(raw)
                values = nan(1, numel(raw));
                for i = 1:numel(raw)
                    if ~isempty(raw{i}) && isnumeric(raw{i})
                        values(i) = raw{i};
                    end
                end
            else
                values = double(raw(:))';
            end
            exact = SAGE.asCellstr(r.exact);
        end

        function [jacobian, latexForm, equilibria] = fluidODEs(rhs, vars, want, url, timeout)
            % FLUIDODES  Jacobian, LaTeX form and equilibria of a fluid drift.
            if nargin < 3 || isempty(want)
                want = {'jacobian', 'latex'};
            end
            if nargin < 4 || isempty(url)
                url = SAGE.require();
            end
            if nargin < 5
                timeout = SAGE.DEFAULT_TIMEOUT;
            end
            payload = struct('rhs', {SAGE.toExpressionList(rhs)}, ...
                'vars', {SAGE.toExpressionList(vars)}, ...
                'want', {want(:)'}, 'timeout_s', timeout);
            r = SAGE.post(url, '/api/v1/fluid/odes', payload, timeout);
            jacobian = {};
            latexForm = {};
            equilibria = {};
            if isfield(r, 'jacobian')
                jacobian = SAGE.asCellMatrix(r.jacobian);
            end
            if isfield(r, 'latex')
                latexForm = SAGE.asCellstr(r.latex);
            end
            if isfield(r, 'equilibria')
                equilibria = r.equilibria;
            end
        end

        % ---------------------------------------------------------------
        % Plumbing
        % ---------------------------------------------------------------

        function url = require()
            % REQUIRE  Resolve a backend or error with how to get one.
            url = SAGE.resolve('auto');
            if isempty(url)
                line_error(mfilename, sprintf(['No symbolic backend is available. Start one with\n', ...
                    '  docker run -d -p 8080:8080 %s\n', ...
                    'point the LINE_SAGE_URL environment variable at a running service, or\n', ...
                    'set options.config.symbolic to its URL.'], SAGE.DOCKER_IMAGES{1}));
            end
        end

        function r = post(url, endpoint, payload, timeout)
            % POST  One JSON request, with the service's own errors raised here.
            if nargin < 4
                timeout = SAGE.DEFAULT_TIMEOUT;
            end
            opts = weboptions('MediaType', 'application/json', ...
                'ContentType', 'json', 'Timeout', max(timeout + 30, 60));
            r = webwrite([SAGE.trimUrl(url), endpoint], payload, opts);
            if ~isfield(r, 'status') || ~strcmp(r.status, 'ok')
                code = 'error';
                msg = 'unspecified error';
                if isfield(r, 'code')
                    code = r.code;
                end
                if isfield(r, 'message')
                    msg = r.message;
                end
                line_error(mfilename, sprintf('line-sage-rest %s failed [%s]: %s', ...
                    endpoint, code, msg));
            end
        end

        function s = trimUrl(url)
            s = regexprep(strtrim(url), '/+$', '');
        end

        function port = freePort()
            % FREEPORT  An ephemeral port nothing is listening on.
            %
            % There is a race between finding the port and docker binding it;
            % it is the same race the JAR and Python clients run, and losing
            % it fails the start loudly rather than silently sharing a port.
            port = 0;
            for attempt = 1:20
                candidate = 20000 + randi(20000);
                [st, ~] = unix(sprintf(...
                    'command -v ss >/dev/null 2>&1 && ss -ltn 2>/dev/null | grep -q ":%d " && echo used', ...
                    candidate));
                if st ~= 0
                    port = candidate;
                    return
                end
            end
            if port == 0
                port = 20000 + randi(20000);
            end
        end

        function [Qcell, symbols] = toExpressionMatrix(Q, symbols)
            % TOEXPRESSIONMATRIX  Generator as a cell array of expressions.
            if iscell(Q)
                Qcell = cellfun(@SAGE.exprToChar, Q, 'UniformOutput', false);
            elseif isnumeric(Q)
                Qcell = arrayfun(@(v) SAGE.numToChar(v), Q, 'UniformOutput', false);
            else
                % sym matrix: char() each entry, keeping ^ which the service
                % reads as exponentiation.
                n = size(Q, 1);
                m = size(Q, 2);
                Qcell = cell(n, m);
                for i = 1:n
                    for j = 1:m
                        Qcell{i, j} = char(Q(i, j));
                    end
                end
                if nargin < 2 || isempty(symbols)
                    symbols = arrayfun(@char, symvar(Q), 'UniformOutput', false);
                end
            end
            % jsonencode maps a cell matrix to a flat array, so hand it a cell
            % of row cells, which becomes the nested array the service wants.
            rows = cell(1, size(Qcell, 1));
            for i = 1:size(Qcell, 1)
                rows{i} = Qcell(i, :);
            end
            Qcell = rows;
            if nargin < 2 || isempty(symbols)
                symbols = {};
            end
            symbols = SAGE.asCellstr(symbols);
        end

        function out = toMeasureStruct(x)
            % TOMEASURESTRUCT  Measures as a struct array with a FIXED field set.
            %
            % jsondecode returns a struct array only when every object in the
            % JSON array carries the SAME fields, and a cell array of structs
            % otherwise. The measure objects are heterogeneous by design: a
            % defined one carries num and den, an undefined one carries reason
            % instead. So the return type flipped between struct array and cell
            % depending on whether any ratio happened to be undefined, which is
            % a property of the MODEL rather than of the call. Normalizing here
            % gives the caller one shape to write against, with '' for whatever
            % the service did not send.
            fields = {'name', 'expr', 'num', 'den', 'reason'};
            if isempty(x)
                out = repmat(cell2struct(repmat({''}, numel(fields), 1), fields, 1), 0, 1);
                return
            end
            if isstruct(x)
                items = num2cell(x(:));
            else
                items = x(:);
            end
            out = repmat(cell2struct(repmat({''}, numel(fields), 1), fields, 1), ...
                numel(items), 1);
            for i = 1:numel(items)
                it = items{i};
                for f = 1:numel(fields)
                    if isstruct(it) && isfield(it, fields{f}) && ~isempty(it.(fields{f}))
                        out(i).(fields{f}) = it.(fields{f});
                    end
                end
            end
        end

        function rows = toExpressionRows(W)
            % TOEXPRESSIONROWS  A possibly NON-SQUARE block as nested cells.
            %
            % toExpressionMatrix is for the generator and assumes a square
            % matrix; a weight block is (measures x states) and rarely square.
            % jsonencode flattens a cell matrix, so the rows are handed over as
            % a cell of row cells, which is the nested array the service wants.
            if isa(W, 'sym')
                C = arrayfun(@char, W, 'UniformOutput', false);
            elseif isnumeric(W)
                C = arrayfun(@(v) SAGE.numToChar(v), W, 'UniformOutput', false);
            else
                C = cellfun(@SAGE.exprToChar, W, 'UniformOutput', false);
            end
            rows = cell(1, size(C, 1));
            for i = 1:size(C, 1)
                rows{i} = C(i, :);
            end
        end

        function out = toExpressionList(exprs)
            % TOEXPRESSIONLIST  Expressions as a cellstr, whatever came in.
            if isa(exprs, 'sym')
                out = arrayfun(@char, exprs(:), 'UniformOutput', false)';
            elseif isnumeric(exprs)
                out = arrayfun(@(v) SAGE.numToChar(v), exprs(:), 'UniformOutput', false)';
            elseif ischar(exprs)
                out = {exprs};
            else
                out = cellfun(@SAGE.exprToChar, exprs(:), 'UniformOutput', false)';
            end
        end

        function s = exprToChar(e)
            if ischar(e)
                s = e;
            elseif isnumeric(e)
                s = SAGE.numToChar(e);
            else
                s = char(e);
            end
        end

        function s = numToChar(v)
            % NUMTOCHAR  Decimal text the service reads as an exact rational.
            %
            % THE TEXT IS THE VALUE, not a rendering of it: the service reads a
            % decimal literal off its source text and converts it to an exact
            % rational, so 0.1 is 1/10 while 0.10000000000000001 is a different
            % number. This prints the SHORTEST decimal that round-trips to the
            % same double, which is what the JAR's Double.toString and the C++
            % decimal_string emit. A fixed '%.17g' does not: on 1/10 it sends
            % 0.10000000000000001 and the same model then yields a different
            % normal form from MATLAB than from the other three codebases.
            if v == round(v) && abs(v) < 1e15
                s = sprintf('%d', round(v));
                return
            end
            for p = 1:17
                s = sprintf('%.*g', p, v);
                if str2double(s) == v
                    return
                end
            end
            s = sprintf('%.17g', v);
        end

        function out = asCellstr(x)
            if isempty(x)
                out = {};
            elseif ischar(x)
                out = {x};
            elseif iscell(x)
                out = reshape(x, 1, []);
            else
                out = num2cell(reshape(x, 1, []));
            end
        end

        function out = asCellMatrix(x)
            % ASCELLMATRIX  Nested JSON array to a cell matrix.
            if iscell(x)
                out = cell(numel(x), 0);
                for i = 1:numel(x)
                    row = SAGE.asCellstr(x{i});
                    out(i, 1:numel(row)) = row;
                end
            else
                out = x;
            end
        end

        function out = toSymOrChar(items)
            % TOSYMORCHAR  Expressions as sym when the toolbox is present.
            %
            % Returning sym where possible keeps every existing caller working
            % unchanged; without the toolbox the same result is still usable
            % as text, which is the whole point of the Sage backend.
            items = SAGE.asCellstr(items);
            if SAGE.hasSymbolicToolbox()
                out = sym(zeros(1, numel(items)));
                for i = 1:numel(items)
                    out(i) = str2sym(items{i});
                end
            else
                out = items;
            end
        end

        function out = scalarSymOrChar(item)
            out = SAGE.toSymOrChar({item});
            if iscell(out)
                out = out{1};
            end
        end

        function bool = hasSymbolicToolbox()
            % HASSYMBOLICTOOLBOX  True if sym objects can be constructed here.
            persistent cached
            if isempty(cached)
                cached = ~isempty(ver('symbolic')) && exist('str2sym', 'file') > 0;
            end
            bool = cached;
        end
    end
end
