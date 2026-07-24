classdef SAGE
    % SAGE  MATLAB-to-SageMath bridge for symbolic analysis.
    %
    % Static helpers that talk to the line-sage-rest service, the same JSON
    % protocol the JAR (jline.api.sym) and native Python
    % (line_solver.api.sym) use. The service is SageMath in a container; see
    % sage/server.py for the endpoints.
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
                if ~SAGE.isReachable(url)
                    url = '';
                end
                return
            end

            env = getenv('LINE_SAGE_URL');
            if ~isempty(env) && SAGE.isReachable(env)
                url = env;
                return
            end

            if ~isempty(cachedUrl) && SAGE.isReachable(cachedUrl)
                url = cachedUrl;
                return
            end

            for p = SAGE.PROBE_PORTS
                candidate = sprintf('http://localhost:%d', p);
                if SAGE.isSageService(candidate)
                    url = candidate;
                    cachedUrl = url;
                    return
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
            if isempty(image)
                url = '';
                return
            end
            url = SAGE.startContainer(image);
            cachedUrl = url;
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

        function url = startContainer(image)
            % STARTCONTAINER  Run the service and wait for it to answer.
            %
            % The container is bound to a free host port, so several MATLAB
            % sessions, or a session next to a hand-started service, do not
            % collide. It is left running: MATLAB has no reliable exit hook,
            % and a warm container is what makes repeated symbolic calls
            % cheap. Stop it with SAGE.stopContainer.
            url = '';
            port = SAGE.freePort();
            name = sprintf('line-sage-rest-%d', port);
            cmd = sprintf('docker run -d --rm --name %s -p %d:8080 %s', name, port, image);
            [st, out] = unix([cmd, ' 2>&1']);
            if st ~= 0
                line_warning(mfilename, 'Could not start %s: %s\n', image, strtrim(out));
                return
            end
            candidate = sprintf('http://localhost:%d', port);
            deadline = tic;
            while toc(deadline) < SAGE.STARTUP_TIMEOUT
                if SAGE.isReachable(candidate)
                    url = candidate;
                    return
                end
                pause(0.5);
            end
            unix(sprintf('docker stop -t 1 %s >/dev/null 2>&1', name));
            line_warning(mfilename, ...
                'Container %s did not answer within %d s.\n', name, SAGE.STARTUP_TIMEOUT);
        end

        function stopContainer()
            % STOPCONTAINER  Stop every container this machine started here.
            if ispc
                return
            end
            unix(['for c in $(docker ps -q --filter name=line-sage-rest-); do ', ...
                'docker stop -t 1 $c >/dev/null 2>&1; done']);
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
            if v == round(v) && abs(v) < 1e15
                s = sprintf('%d', round(v));
            else
                s = sprintf('%.17g', v);
            end
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
