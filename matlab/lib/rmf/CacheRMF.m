% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

classdef CacheRMF
    % CacheRMF Multi-list cache with RANDOM(m) replacement as a DDPP.
    %
    % Models a cache with h lists of capacities m = [m_1, ..., m_h] and n
    % items with popularity distribution p = [p_1, ..., p_n]. The state
    % space tracks in which list each item resides (or if it is outside the
    % cache).
    %
    % State vector layout (1-based MATLAB indexing):
    %   x(i + k * n) = density of item i in list k
    %   where i = 1..n, k = 0..h (k=0 means outside cache, k=1..h are
    %   cache lists)
    %
    % Transition dynamics (RANDOM(m) replacement):
    %   When item i is requested (rate p_i) and item i is in list k:
    %   - A uniformly random item j from list k+1 is displaced to list k
    %   - Item i moves from list k to list k+1
    %
    % Reference:
    %   N. Gast, "Expected Values Estimated via Mean-Field Approximation
    %   are 1/N-Accurate", Proc. ACM Meas. Anal. Comput. Syst., 2017.

    properties
        p               % double array (1 x n), item popularity distribution
        m               % int array (1 x h), list capacities
        numberOfItems   % int, number of items (n)
        numberOfLists   % int, number of cache lists (h)
        modelDimension  % int, total state dimension n*(h+1)
        x0              % double array (1 x modelDimension), initial state
    end

    methods
        function obj = CacheRMF(p, m)
            % CacheRMF Constructor.
            %
            % Parameters:
            %   p - double array (1 x n), item request probabilities
            %   m - int array (1 x h), capacity of each cache list
            %
            % Builds initial state: first m(1) items in list 1, next m(2)
            % in list 2, etc. Remaining items outside cache (list 0).

            obj.p = p(:)';  % ensure row vector
            obj.m = m(:)';  % ensure row vector
            obj.numberOfItems = length(p);
            obj.numberOfLists = length(m);
            obj.modelDimension = obj.numberOfItems * (obj.numberOfLists + 1);

            n = obj.numberOfItems;
            h = obj.numberOfLists;

            % Build initial state
            obj.x0 = zeros(1, obj.modelDimension);
            objIdx = 1;  % 1-based item counter
            for k = 1:h
                for c = 1:m(k)
                    if objIdx <= n
                        obj.x0(obj.index(objIdx, k)) = 1.0;
                        objIdx = objIdx + 1;
                    end
                end
            end
            for i = objIdx:n
                obj.x0(obj.index(i, 0)) = 1.0;
            end
        end

        function idx = index(obj, i, k)
            % index Map (item i, list k) to flat state vector index.
            %
            % Parameters:
            %   i - item index (1-based, 1..n), scalar or vector
            %   k - list index (0-based, 0..h), scalar
            %
            % Returns:
            %   idx - 1-based flat index into state vector

            idx = i + k * obj.numberOfItems;
        end

        function hr = hitRate(obj, x, listNumber)
            % hitRate Compute hit rate contribution from a specific list.
            %
            % Parameters:
            %   x          - state vector (1 x modelDimension) or column vector
            %   listNumber - list index (0 = outside cache, 1..h = cache lists)
            %
            % Returns:
            %   hr - sum of p(i) * x(index(i, listNumber)) over all items

            indices = obj.index(1:obj.numberOfItems, listNumber);
            hr = dot(obj.p(:), x(indices(:)));
        end

        function dX = drift(obj, x)
            % drift Compute the mean field drift F(x).
            %
            % The drift for RANDOM(m) replacement is:
            %   dx(i,k)/dt = -p(i)*x(i,k) + hitRate(k)*x(i,k+1)/m(k)
            %   dx(i,k+1)/dt = p(i)*x(i,k) - hitRate(k)*x(i,k+1)/m(k)
            %
            % Parameters:
            %   x - state vector (modelDimension x 1) or (1 x modelDimension)
            %
            % Returns:
            %   dX - drift vector, same shape as x

            x = x(:)';  % ensure row vector for internal computation
            n = obj.numberOfItems;
            h = obj.numberOfLists;

            hitRates = zeros(1, h + 1);
            for k = 0:h
                hitRates(k + 1) = obj.hitRate(x, k);
            end

            dX = zeros(1, obj.modelDimension);
            for i = 1:n
                for k = 0:(h - 1)
                    ik = obj.index(i, k);
                    ik1 = obj.index(i, k + 1);
                    flow = obj.p(i) * x(ik) - hitRates(k + 1) * x(ik1) / obj.m(k + 1);
                    dX(ik) = dX(ik) - flow;
                    dX(ik1) = dX(ik1) + flow;
                end
            end
            dX = dX(:);  % return column vector (for ODE solvers)
        end

        function Fp = jacobian(obj, x)
            % jacobian Compute Jacobian dF/dx at state x.
            %
            % Parameters:
            %   x - state vector (1 x modelDimension)
            %
            % Returns:
            %   Fp - Jacobian matrix (modelDimension x modelDimension)

            x = x(:)';
            n = obj.numberOfItems;
            h = obj.numberOfLists;
            dim = obj.modelDimension;

            hitRates = zeros(1, h + 1);
            for k = 0:h
                hitRates(k + 1) = obj.hitRate(x, k);
            end

            Fp = zeros(dim, dim);
            for i = 1:n
                for k = 0:(h - 1)
                    ik = obj.index(i, k);
                    ik1 = obj.index(i, k + 1);
                    mk = obj.m(k + 1);  % m is 1-indexed, k is 0-based list

                    % Direct rate terms
                    Fp(ik, ik) = Fp(ik, ik) - obj.p(i);
                    Fp(ik1, ik) = Fp(ik1, ik) + obj.p(i);
                    Fp(ik, ik1) = Fp(ik, ik1) + hitRates(k + 1) / mk;
                    Fp(ik1, ik1) = Fp(ik1, ik1) - hitRates(k + 1) / mk;

                    % Indirect terms via hit rate dependence
                    for j = 1:n
                        jk = obj.index(j, k);
                        jk1 = obj.index(j, k + 1);
                        % d(hitRate[k])/d(x[j,k+1]) = p[j] for list k+1
                        % but here hitRates(k+1) is hitRate for list k
                        Fp(ik, jk1) = Fp(ik, jk1) - obj.p(i) * x(ik) / mk;
                        Fp(ik1, jk1) = Fp(ik1, jk1) + obj.p(i) * x(ik) / mk;
                        Fp(ik, jk) = Fp(ik, jk) + obj.p(j) * x(ik1) / mk;
                        Fp(ik1, jk) = Fp(ik1, jk) - obj.p(j) * x(ik1) / mk;
                    end
                end
            end
        end

        function Fpp = hessian(obj, x) %#ok<INUSD>
            % hessian Compute Hessian d^2F/dx^2.
            %
            % The Hessian is constant (drift is quadratic in x), so x is
            % unused but kept for interface consistency.
            %
            % Parameters:
            %   x - state vector (unused)
            %
            % Returns:
            %   Fpp - Hessian tensor (modelDimension x modelDimension x modelDimension)

            n = obj.numberOfItems;
            h = obj.numberOfLists;
            dim = obj.modelDimension;
            Fpp = zeros(dim, dim, dim);

            for i = 1:n
                for k = 0:(h - 1)
                    ik = obj.index(i, k);
                    ik1 = obj.index(i, k + 1);
                    mk = obj.m(k + 1);
                    for j = 1:n
                        if j ~= i
                            jk = obj.index(j, k);
                            jk1 = obj.index(j, k + 1);
                            % d^2 F[ik] / (d x[jk] d x[ik1]) = p[j]/mk
                            Fpp(ik, jk, ik1) = Fpp(ik, jk, ik1) + obj.p(j) / mk;
                            Fpp(ik, ik1, jk) = Fpp(ik, ik1, jk) + obj.p(j) / mk;
                            % d^2 F[ik] / (d x[jk1] d x[ik]) = -p[i]/mk
                            Fpp(ik, jk1, ik) = Fpp(ik, jk1, ik) - obj.p(i) / mk;
                            Fpp(ik, ik, jk1) = Fpp(ik, ik, jk1) - obj.p(i) / mk;
                            % Symmetric for ik1
                            Fpp(ik1, jk, ik1) = Fpp(ik1, jk, ik1) - obj.p(j) / mk;
                            Fpp(ik1, ik1, jk) = Fpp(ik1, ik1, jk) - obj.p(j) / mk;
                            Fpp(ik1, jk1, ik) = Fpp(ik1, jk1, ik) + obj.p(i) / mk;
                            Fpp(ik1, ik, jk1) = Fpp(ik1, ik, jk1) + obj.p(i) / mk;
                        end
                    end
                end
            end
        end

        function Q = noiseMatrix(obj, x)
            % noiseMatrix Compute noise intensity matrix Q(x) for the DDPP.
            %
            % Q(a,b) = sum_ell ell(a) * ell(b) * beta_ell(x)
            % where each transition ell is a swap between items i and j
            % across lists k and k+1, with rate p(i)*x(i,k)*x(j,k+1)/m(k).
            %
            % Parameters:
            %   x - state vector (1 x modelDimension)
            %
            % Returns:
            %   Q - noise matrix (modelDimension x modelDimension)

            x = x(:)';
            n = obj.numberOfItems;
            h = obj.numberOfLists;
            dim = obj.modelDimension;
            Q = zeros(dim, dim);

            signs = [-1, 1, 1, -1];
            for i = 1:n
                for k = 0:(h - 1)
                    mk = obj.m(k + 1);
                    ik = obj.index(i, k);
                    ik1 = obj.index(i, k + 1);
                    for j = 1:n
                        jk = obj.index(j, k);
                        jk1 = obj.index(j, k + 1);
                        rate = obj.p(i) * x(ik) * x(jk1) / mk;
                        indices = [ik, jk, ik1, jk1];
                        for ia = 1:4
                            for ib = 1:4
                                Q(indices(ia), indices(ib)) = Q(indices(ia), indices(ib)) + rate * signs(ia) * signs(ib);
                            end
                        end
                    end
                end
            end
        end

        function pi = fixedPoint(obj, tmax)
            % fixedPoint Compute mean field fixed point by ODE integration.
            %
            % Integrates dx/dt = F(x) until steady state using ode15s.
            %
            % Parameters:
            %   tmax - maximum integration time (default: 10000)
            %
            % Returns:
            %   pi - fixed point state vector (1 x modelDimension)

            if nargin < 2
                tmax = 10000;
            end
            opts = odeset('RelTol', 1e-8, 'AbsTol', 1e-10);
            [~, Y] = ode15s(@(t, x) obj.drift(x), [0, tmax], obj.x0', opts);
            %[~, Y] = lsoda_solve(@(t, x) obj.drift(x), [0, tmax], obj.x0', opts);
            pi = Y(end, :);
        end

        function [P, Pinv, r] = dimensionReduction(obj, Fp)
            % dimensionReduction SVD-based dimension reduction for singular Jacobian.
            %
            % The Jacobian is singular because item populations are conserved
            % (sum over lists for each item = 1). This method finds a change
            % of basis that separates the rank-deficient directions.
            %
            % Parameters:
            %   Fp - Jacobian matrix (modelDimension x modelDimension)
            %
            % Returns:
            %   P    - change-of-basis matrix
            %   Pinv - inverse of P
            %   r    - rank of Fp

            dim = obj.modelDimension;
            n = obj.numberOfItems;
            h = obj.numberOfLists;

            r = rank(Fp);

            % Build change-of-basis: first r rows are independent coordinates,
            % remaining rows span the null space of Fp
            C = zeros(dim, dim);
            d = 1;
            for lIdx = 0:h
                for i = 1:(n - 1)
                    C(d, obj.index(i, lIdx)) = 1.0;
                    d = d + 1;
                end
            end

            [U, ~, ~] = svd(Fp);
            C((r + 1):end, :) = U(:, (r + 1):end)';
            Pinv = inv(C); %#ok<MINV>
            P = C;
        end

        function [Fp_r, Fpp_r, Q_r, P, Pinv, r] = reduceFpFppQ(obj, Fp, Fpp, Q)
            % reduceFpFppQ Apply dimension reduction to Fp, Fpp, Q.
            %
            % Projects the Jacobian, Hessian, and noise matrix onto the
            % non-singular subspace identified by dimension reduction.
            %
            % Parameters:
            %   Fp  - Jacobian matrix (dim x dim)
            %   Fpp - Hessian tensor (dim x dim x dim)
            %   Q   - noise matrix (dim x dim)
            %
            % Returns:
            %   Fp_r  - reduced Jacobian (r x r)
            %   Fpp_r - reduced Hessian (r x r x r)
            %   Q_r   - reduced noise matrix (r x r)
            %   P     - change-of-basis matrix
            %   Pinv  - inverse of P
            %   r     - rank

            [P, Pinv, r] = obj.dimensionReduction(Fp);

            % Reduced Jacobian
            Fp_full = P * Fp * Pinv;
            Fp_r = Fp_full(1:r, 1:r);

            % Reduced Hessian via tensor contraction
            % Fpp_r(a,b,c) = sum_{i,j,k} P(a,i) * Fpp(i,j,k) * Pinv(j,b) * Pinv(k,c)
            dim = obj.modelDimension;
            Fpp_r = zeros(r, r, r);
            % First contraction: T1(a,j,k) = sum_i P(a,i) * Fpp(i,j,k)
            T1 = zeros(r, dim, dim);
            for a = 1:r
                for j = 1:dim
                    for k = 1:dim
                        val = 0;
                        for i = 1:dim
                            val = val + P(a, i) * Fpp(i, j, k);
                        end
                        T1(a, j, k) = val;
                    end
                end
            end
            % Second contraction: T2(a,b,k) = sum_j T1(a,j,k) * Pinv(j,b)
            T2 = zeros(r, r, dim);
            for a = 1:r
                for b = 1:r
                    for k = 1:dim
                        val = 0;
                        for j = 1:dim
                            val = val + T1(a, j, k) * Pinv(j, b);
                        end
                        T2(a, b, k) = val;
                    end
                end
            end
            % Third contraction: Fpp_r(a,b,c) = sum_k T2(a,b,k) * Pinv(k,c)
            for a = 1:r
                for b = 1:r
                    for c = 1:r
                        val = 0;
                        for k = 1:dim
                            val = val + T2(a, b, k) * Pinv(k, c);
                        end
                        Fpp_r(a, b, c) = val;
                    end
                end
            end

            % Reduced noise matrix
            Q_full = P * Q * P';
            Q_r = Q_full(1:r, 1:r);
        end

        function [V, W] = expandVW(obj, V_r, W_r, Pinv, r) %#ok<INUSL>
            % expandVW Expand reduced V, W back to full dimension.
            %
            % Parameters:
            %   V_r  - reduced correction vector (r x 1)
            %   W_r  - reduced covariance matrix (r x r)
            %   Pinv - inverse change-of-basis matrix
            %   r    - rank
            %
            % Returns:
            %   V - full correction vector (modelDimension x 1)
            %   W - full covariance matrix (modelDimension x modelDimension)

            V = Pinv(:, 1:r) * V_r;
            W = Pinv(:, 1:r) * W_r * Pinv(:, 1:r)';
        end

        function [pi, V, VW] = meanFieldExpansionSteadyState(obj, order)
            % meanFieldExpansionSteadyState Compute refined mean field steady-state.
            %
            % Computes the mean field fixed point pi and the 1/N correction V
            % using the Lyapunov equation approach with dimension reduction.
            %
            % The refined approximation for a system of N items is:
            %   E[X] ~ pi + V/N + O(1/N^2)
            %
            % Parameters:
            %   order - expansion order (0 = plain MF, 1 = with 1/N correction)
            %           default: 1
            %
            % Returns:
            %   pi - mean field fixed point (1 x modelDimension)
            %   V  - first-order correction (modelDimension x 1)
            %   VW - cell {V, W} with correction and covariance

            if nargin < 2
                order = 1;
            end

            pi = obj.fixedPoint();

            if order == 0
                V = zeros(obj.modelDimension, 1);
                W = zeros(obj.modelDimension);
                VW = {V, W};
                return;
            end

            Fp = obj.jacobian(pi);
            Fpp = obj.hessian(pi);
            Q = obj.noiseMatrix(pi);

            % Dimension reduction: project onto non-singular subspace
            [Fp_r, Fpp_r, Q_r, ~, Pinv, r] = obj.reduceFpFppQ(Fp, Fpp, Q);

            % Solve Lyapunov equation in reduced space:
            %   Fp_r * W_r + W_r * Fp_r' + Q_r = 0
            % MATLAB lyap(A, Q) solves A*X + X*A' + Q = 0
            W_r = lyap(Fp_r, Q_r);

            % First-order correction: Hessian contraction
            % C_r(a) = sum_{b,c} Fpp_r(a,b,c) * W_r(b,c)
            C_r = zeros(r, 1);
            for a = 1:r
                C_r(a) = sum(sum(squeeze(Fpp_r(a, :, :)) .* W_r));
            end
            V_r = -Fp_r \ (C_r / 2);

            % Expand back to full dimension
            [V, W] = obj.expandVW(V_r, W_r, Pinv, r);
            VW = {V, W};
        end

        function [T, X, Vt, Wt] = meanFieldExpansionTransient(obj, time, n_points, order)
            % meanFieldExpansionTransient Compute refined mean field transient expansion.
            %
            % Integrates the coupled ODE system for (X, V, W) where:
            %   X(t): mean field trajectory
            %   V(t): 1/N correction trajectory
            %   W(t): covariance trajectory
            %
            % Parameters:
            %   time     - maximum integration time (default: 50)
            %   n_points - number of output time points (default: 200)
            %   order    - expansion order (0 or 1, default: 1)
            %
            % Returns:
            %   T  - time points (n_points x 1)
            %   X  - mean field trajectory (n_points x modelDimension)
            %   Vt - correction trajectory (n_points x modelDimension)
            %   Wt - covariance trajectory (n_points x modelDimension x modelDimension)

            if nargin < 2 || isempty(time), time = 50; end
            if nargin < 3 || isempty(n_points), n_points = 200; end
            if nargin < 4 || isempty(order), order = 1; end

            dim = obj.modelDimension;
            T = linspace(0, time, n_points)';

            if order == 0
                odeopt = odeset('RelTol', 1e-6, 'AbsTol', 1e-10);
                [~, X] = ode15s(@(t, x) obj.drift(x), T, obj.x0', odeopt);
                %[~, X] = lsoda_solve(@(t, x) obj.drift(x), T, obj.x0', odeopt);
                Vt = zeros(n_points, dim);
                Wt = zeros(n_points, dim, dim);
                return;
            end

            % Coupled ODE: state = [X (dim), V (dim), W_flat (dim^2)]
            total_dim = dim + dim + dim * dim;
            y0 = zeros(total_dim, 1);
            y0(1:dim) = obj.x0(:);

            function dy = coupled_rhs(~, y)
                x = y(1:dim);
                v = y((dim+1):(2*dim));
                w = reshape(y((2*dim+1):end), dim, dim);

                F = obj.drift(x);
                Fp = obj.jacobian(x(:)');
                Fpp = obj.hessian(x(:)');
                Qmat = obj.noiseMatrix(x(:)');

                dx = F(:);
                dv = Fp * v;
                % Add Hessian contraction: 0.5 * tensordot(Fpp, w, axes=([1,2],[0,1]))
                for a = 1:dim
                    dv(a) = dv(a) + 0.5 * sum(sum(squeeze(Fpp(a, :, :)) .* w));
                end
                dw = Fp * w + w * Fp' + Qmat;

                dy = [dx; dv; dw(:)];
            end

            odeopt = odeset('RelTol', 1e-6, 'AbsTol', 1e-10);
            [~, Y] = ode15s(@coupled_rhs, T, y0, odeopt);
            %[~, Y] = lsoda_solve(@coupled_rhs, T, y0, odeopt);

            X = Y(:, 1:dim);
            Vt = Y(:, (dim+1):(2*dim));
            W_flat = Y(:, (2*dim+1):end);
            Wt = reshape(W_flat', dim, dim, []);
            Wt = permute(Wt, [3, 1, 2]);
        end

        function hr = hitRatesAll(obj, x)
            % hitRatesAll Compute hit rates for all lists.
            %
            % Parameters:
            %   x - state vector (1 x modelDimension)
            %
            % Returns:
            %   hr - array (1 x numberOfLists+1) with hit rate per list

            h = obj.numberOfLists;
            hr = zeros(1, h + 1);
            for k = 0:h
                hr(k + 1) = obj.hitRate(x, k);
            end
        end

        function out = convertTo2D(obj, x)
            % convertTo2D Convert flat state vector to (n x h+1) matrix.
            %
            % Parameters:
            %   x - state vector (1 x modelDimension) or matrix
            %       (T x modelDimension) for trajectories
            %
            % Returns:
            %   out - array (n x h+1) or (T x n x h+1)

            n = obj.numberOfItems;
            h = obj.numberOfLists;
            if isvector(x)
                x = x(:)';
                out = zeros(n, h + 1);
                for i = 1:n
                    for k = 0:h
                        out(i, k + 1) = x(obj.index(i, k));
                    end
                end
            else
                T = size(x, 1);
                out = zeros(T, n, h + 1);
                for t = 1:T
                    for i = 1:n
                        for k = 0:h
                            out(t, i, k + 1) = x(t, obj.index(i, k));
                        end
                    end
                end
            end
        end
    end
end
