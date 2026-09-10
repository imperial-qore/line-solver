classdef CME < ME
    % Concentrated Matrix Exponential (CME) distribution
    %
    % A CME is the matrix-exponential distribution of odd order 2*n+1 whose
    % squared coefficient of variation is (numerically) minimal for that order,
    % from the tables of Horvath, Horvath and Telek. Its SCV decays as O(1/n^2)
    % and therefore goes far below the Erlang bound 1/order attainable by a
    % phase-type distribution of the same order: order 101 gives SCV 3.9e-4,
    % where Erlang-101 gives 9.9e-3.
    %
    % The density of the unit-mean CME with n harmonic terms is
    %
    %   f(x) = mu1*exp(-mu1*x)*(c + sum_k [a_k*cos(k*w*mu1*x) + b_k*sin(k*w*mu1*x)])
    %
    % with w = omega, which is exactly alpha*expm(A*x)*(-A*e) for the
    % block-diagonal A = blkdiag(-mu1, mu1*[-1 -k*w; k*w -1], k=1..n). The
    % parameters a, b, c, omega and mu1 are read from the same iltcme.json
    % table used by the CME inverse Laplace transform (matlab_ilt).
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        cmeMean;   % Mean of the distribution
        cmeOrder;  % Number of phases, an odd integer 2*n+1
    end

    methods
        function self = CME(mean, order)
            % SELF = CME(MEAN, ORDER)
            % Create a concentrated matrix-exponential distribution
            %
            % @param mean Mean of the distribution (positive)
            % @param order Number of phases, an odd integer 2*n+1 with n in the table
            % @return self CME distribution instance

            [alpha, A] = CME.representation(order);

            if ~isscalar(mean) || ~isfinite(mean) || mean <= 0
                line_error(mfilename, 'CME mean must be a positive finite number.');
            end

            % The CME density is nonnegative by construction, so the density
            % scan of ME is skipped: it can never fire and it costs O(1e5)
            % propagations of a (2n+1)-square matrix, prohibitive at high order.
            self@ME(alpha, A/mean, false);

            self.cmeMean = mean;
            self.cmeOrder = double(order);

            % The process type stays ME: a CME is a matrix-exponential
            % representation, so every solver gate, sn.procid entry and JSON key
            % that accepts ME accepts it unchanged.
            self.obj = jline.lang.processes.CME(mean, double(order));

            % see _kb/04-networkstruct.md (CME.m representation invariants)
        end

        function o = getOrder(self)
            % O = GETORDER()
            % Get the CME order, i.e. the number of phases
            o = self.cmeOrder;
        end
    end

    methods(Static)
        function params = table()
            % PARAMS = TABLE()
            % Load and cache the CME parameter table from iltcme.json
            %
            % The same table backs matlab_ilt, which caches it in the global
            % cmeParams, so the two share one decode per session.

            global cmeParams;
            if isempty(cmeParams)
                cmeParams = jsondecode(fileread('iltcme.json'));
            end
            params = cmeParams;
        end

        function entry = tableEntry(order)
            % ENTRY = TABLEENTRY(ORDER)
            % Select the CME table entry realizing the given number of phases
            %
            % The table is keyed by the number of harmonic terms n, so an order
            % of 2*n+1 phases maps to the entries with that n. Several entries
            % can share an n (the 'full' and 'approx' optimizations), and the
            % most concentrated one is taken, matching the selection rule of the
            % CME inverse Laplace transform.

            order = double(order);
            if ~isscalar(order) || order < 3 || mod(order,2) == 0 || order ~= fix(order)
                line_error(mfilename, sprintf('CME order must be an odd integer of the form 2*n+1 with n >= 1, got %g.', order));
            end
            n = (order-1)/2;

            params = CME.table();
            entry = [];
            bestcv2 = Inf;
            for i = 1:numel(params)
                if iscell(params)
                    cand = params{i};
                else
                    cand = params(i);
                end
                if cand.n == n && cand.cv2 < bestcv2
                    entry = cand;
                    bestcv2 = cand.cv2;
                end
            end

            if isempty(entry)
                orders = CME.getSupportedOrders();
                [~,j] = min(abs(orders-order));
                line_error(mfilename, sprintf('No tabulated CME of order %d; the nearest available order is %d. Use CME.getSupportedOrders() for the full list.', order, orders(j)));
            end
        end

        function orders = getSupportedOrders()
            % ORDERS = GETSUPPORTEDORDERS()
            % Get the sorted list of CME orders (phase counts) in the table

            params = CME.table();
            ns = zeros(numel(params),1);
            for i = 1:numel(params)
                if iscell(params)
                    ns(i) = params{i}.n;
                else
                    ns(i) = params(i).n;
                end
            end
            orders = unique(2*ns+1)';
        end

        function scv = getMinSCV(order)
            % SCV = GETMINSCV(ORDER)
            % Get the tabulated minimal SCV attained by a CME of the given order

            entry = CME.tableEntry(order);
            scv = entry.cv2;
        end

        function cme = fitMeanAndSCV(mean, scv)
            % CME = FITMEANANDSCV(MEAN, SCV)
            % Create the lowest-order CME with the given mean and SCV at most scv
            %
            % @param mean Target mean
            % @param scv Target squared coefficient of variation, an upper bound.
            %        The lowest tabulated order whose minimal SCV does not exceed
            %        it is selected, so the result is at least as concentrated as
            %        requested.
            % @return cme CME of that order rescaled to the requested mean

            params = CME.table();
            bestorder = Inf;
            mincv2 = Inf;
            maxn = 0;
            for i = 1:numel(params)
                if iscell(params)
                    cand = params{i};
                else
                    cand = params(i);
                end
                mincv2 = min(mincv2, cand.cv2);
                maxn = max(maxn, cand.n);
                if cand.cv2 <= scv
                    bestorder = min(bestorder, 2*cand.n+1);
                end
            end

            if ~isfinite(bestorder)
                line_error(mfilename, sprintf('No tabulated CME reaches SCV %g; the most concentrated entry has SCV %g at order %d.', scv, mincv2, 2*maxn+1));
            end

            cme = CME(mean, bestorder);
        end

        function [alpha, A, scv] = representation(order)
            % [ALPHA, A, SCV] = REPRESENTATION(ORDER)
            % Build the unit-mean (alpha, A) matrix-exponential form of a CME
            %
            % A = blkdiag(-mu1, mu1*[-1 -k*w; k*w -1], k=1..n) reproduces the
            % exponential envelope in its first phase and the k-th harmonic in
            % its k-th 2x2 rotation block, since expm(mu1*x*[-1 -kw; kw -1]) is
            % exp(-mu1*x) times the rotation by k*w*mu1*x. The entries of alpha
            % follow by matching -alpha*expm(A*x)*A*e term by term: with
            % wk = k*w and d = 2*(1+wk^2),
            %
            %   alpha(1)     = c
            %   alpha(2k)    = ((1+wk)*a_k - (1-wk)*b_k)/d
            %   alpha(2k+1)  = ((1-wk)*a_k + (1+wk)*b_k)/d
            %
            % The result has unit mean and sums to one, as any (alpha, A) whose
            % density integrates to one must.

            entry = CME.tableEntry(order);
            a = entry.a(:);
            b = entry.b(:);
            c = entry.c;
            mu1 = entry.mu1;
            w = entry.omega;
            n = entry.n;

            sz = 2*n+1;
            A = zeros(sz, sz);
            alpha = zeros(1, sz);
            A(1,1) = -mu1;
            alpha(1) = c;
            for k = 1:n
                i = 2*k;
                wk = k*w;
                A(i,i) = -mu1;
                A(i,i+1) = -wk*mu1;
                A(i+1,i) = wk*mu1;
                A(i+1,i+1) = -mu1;
                d = 2*(1+wk^2);
                alpha(i) = ((1+wk)*a(k) - (1-wk)*b(k))/d;
                alpha(i+1) = ((1-wk)*a(k) + (1+wk)*b(k))/d;
            end

            % see _kb/04-networkstruct.md (CME.m representation invariants)
            alpha = alpha / sum(alpha);

            scv = entry.cv2;
        end
    end
end
