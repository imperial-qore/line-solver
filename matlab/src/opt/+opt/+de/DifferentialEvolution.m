classdef DifferentialEvolution < handle
    % DifferentialEvolution  Self-contained port of scipy's
    % differential_evolution that reproduces its trajectory bit-for-bit for a
    % given integer seed, using opt.de.NumpyRandomState for all draws.
    %
    % Covers the configuration used by line-opt: binomial strategies (default
    % best1bin), updating='immediate', dithered mutation, latinhypercube
    % initialization, polish=false, and penalty-based constraints handled inside
    % the objective. The exact numpy draw order is preserved: LHS init (uniform
    % grid then per-column permutation), a per-generation dither uniform, and per
    % candidate randint (fill point), shuffle (sample selection) and uniform
    % (crossover), plus a data-dependent uniform in the bound-repair step.

    properties
        objective          % function handle: energy = f(x)
        low                % 1xd lower bounds
        high               % 1xd upper bounds
        paramCount
        strategy
        popsizeMult
        maxiter
        ditherLow
        ditherHigh
        recombination
        tol
        atol
        callback           % function handle (bestX, nit) -> logical stop, or []
        rng                % opt.de.NumpyRandomState
        numMembers
        population         % numMembers x d, scaled to [0,1]
        energies           % 1 x numMembers
        randomIndex        % persistent shuffle buffer (0-based values)
        scale
        nfev
        bestPerGen         % cell array of 1xd vectors
    end

    methods
        function obj = DifferentialEvolution(objective, low, high, strategy, ...
                popsizeMult, maxiter, ditherLow, ditherHigh, recombination, tol, seed)
            obj.objective = objective;
            obj.low = low(:).';
            obj.high = high(:).';
            obj.paramCount = numel(low);
            obj.strategy = strategy;
            obj.popsizeMult = popsizeMult;
            obj.maxiter = maxiter;
            obj.ditherLow = min(ditherLow, ditherHigh);
            obj.ditherHigh = max(ditherLow, ditherHigh);
            obj.recombination = recombination;
            obj.tol = tol;
            obj.atol = 0.0;
            obj.callback = [];
            obj.scale = ditherLow;
            obj.rng = opt.de.NumpyRandomState(seed);
            obj.bestPerGen = {};
            obj.initPopulationLhs();
        end

        function s1 = scaleArg1(obj, j)
            s1 = 0.5 * (obj.low(j) + obj.high(j));
        end

        function s2 = scaleArg2(obj, j)
            s2 = abs(obj.low(j) - obj.high(j));
        end

        function out = scaleParameters(obj, trial)
            out = zeros(1, obj.paramCount);
            for j = 1:obj.paramCount
                out(j) = obj.scaleArg1(j) + (trial(j) - 0.5) * obj.scaleArg2(j);
            end
        end

        function initPopulationLhs(obj)
            ebCount = sum(obj.low == obj.high);
            obj.numMembers = max(5, obj.popsizeMult * max(1, obj.paramCount - ebCount));
            segsize = 1.0 / obj.numMembers;
            flat = obj.rng.uniformN(0, 1, obj.numMembers * obj.paramCount);
            samples = zeros(obj.numMembers, obj.paramCount);
            p = 1;
            for i = 1:obj.numMembers
                offset = (i - 1) / obj.numMembers;
                for j = 1:obj.paramCount
                    samples(i, j) = segsize * flat(p) + offset;
                    p = p + 1;
                end
            end
            obj.population = zeros(obj.numMembers, obj.paramCount);
            for j = 1:obj.paramCount
                order = obj.rng.permutation(obj.numMembers);  % 0-based
                for i = 1:obj.numMembers
                    obj.population(i, j) = samples(order(i) + 1, j);
                end
            end
            obj.energies = inf(1, obj.numMembers);
            obj.randomIndex = 0:(obj.numMembers - 1);
            obj.nfev = 0;
        end

        function b = best(obj)
            b = obj.scaleParameters(obj.population(1, :));
        end

        function calculateInitialEnergies(obj)
            for i = 1:obj.numMembers
                obj.energies(i) = obj.objective(obj.scaleParameters(obj.population(i, :)));
                obj.nfev = obj.nfev + 1;
            end
        end

        function promoteLowestEnergy(obj)
            [~, l] = min(obj.energies);
            if l ~= 1
                te = obj.energies(1); obj.energies(1) = obj.energies(l); obj.energies(l) = te;
                tp = obj.population(1, :); obj.population(1, :) = obj.population(l, :); obj.population(l, :) = tp;
            end
        end

        function tf = converged(obj)
            if any(isinf(obj.energies))
                tf = false;
                return;
            end
            tf = std(obj.energies, 1) <= obj.atol + obj.tol * abs(mean(obj.energies));
        end

        function out = selectSamples(obj, candidate, numberSamples)
            % candidate is 0-based member index
            obj.randomIndex = obj.rng.shuffle(obj.randomIndex);
            pick = obj.randomIndex(1:numberSamples + 1);
            out = zeros(1, numberSamples);
            k = 0;
            for i = 1:numel(pick)
                if k >= numberSamples
                    break;
                end
                if pick(i) ~= candidate
                    k = k + 1;
                    out(k) = pick(i);
                end
            end
        end

        function b = bprime(obj, candidate, s)
            % candidate and s are 0-based member indices; convert to 1-based rows
            P = obj.population;
            c = candidate + 1;
            r = s + 1;
            switch obj.strategy
                case {'rand1bin', 'rand1exp'}
                    b = P(r(1), :) + obj.scale * (P(r(2), :) - P(r(3), :));
                case {'randtobest1bin', 'randtobest1exp'}
                    b = P(r(1), :);
                    b = b + obj.scale * (P(1, :) - b);
                    b = b + obj.scale * (P(r(2), :) - P(r(3), :));
                case {'currenttobest1bin', 'currenttobest1exp'}
                    b = P(c, :) + obj.scale * (P(1, :) - P(c, :) + P(r(1), :) - P(r(2), :));
                case {'best2bin', 'best2exp'}
                    b = P(1, :) + obj.scale * (P(r(1), :) + P(r(2), :) - P(r(3), :) - P(r(4), :));
                case {'rand2bin', 'rand2exp'}
                    b = P(r(1), :) + obj.scale * (P(r(2), :) + P(r(3), :) - P(r(4), :) - P(r(5), :));
                otherwise   % best1bin / best1exp
                    b = P(1, :) + obj.scale * (P(r(1), :) - P(r(2), :));
            end
        end

        function trial = mutate(obj, candidate)
            % candidate is 0-based member index
            fillPoint = obj.rng.randint(0, obj.paramCount);   % 0-based
            samples = obj.selectSamples(candidate, 5);
            bp = obj.bprime(candidate, samples);
            trial = obj.population(candidate + 1, :);
            cross = obj.rng.uniformN(0, 1, obj.paramCount);
            doCross = cross < obj.recombination;
            doCross(fillPoint + 1) = true;
            trial(doCross) = bp(doCross);
        end

        function trial = ensureConstraint(obj, trial)
            mask = (trial > 1) | (trial < 0);
            oob = sum(mask);
            if oob > 0
                repl = obj.rng.uniformN(0, 1, oob);
                trial(mask) = repl;
            end
        end

        function result = solve(obj)
            warningFlag = false;
            if any(isinf(obj.energies))
                obj.calculateInitialEnergies();
                obj.promoteLowestEnergy();
            end
            nit = 0;
            for nit = 1:obj.maxiter
                obj.next();
                if ~isempty(obj.callback)
                    stop = obj.callback(obj.best(), nit);
                    if stop
                        warningFlag = true;
                    end
                end
                obj.bestPerGen{end+1} = obj.best(); %#ok<AGROW>
                if warningFlag || obj.converged()
                    break;
                end
            end
            if nit >= obj.maxiter
                warningFlag = warningFlag || (nit == obj.maxiter && ~obj.converged());
            end
            result = struct();
            result.x = obj.best();
            result.fun = obj.energies(1);
            result.nit = nit;
            result.nfev = obj.nfev;
            result.success = ~warningFlag;
            result.bestPerGen = obj.bestPerGen;
        end

        function next(obj)
            obj.scale = obj.rng.uniformScalar(obj.ditherLow, obj.ditherHigh);
            for candidate = 0:(obj.numMembers - 1)
                trial = obj.mutate(candidate);
                trial = obj.ensureConstraint(trial);
                parameters = obj.scaleParameters(trial);
                energy = obj.objective(parameters);
                obj.nfev = obj.nfev + 1;
                if energy <= obj.energies(candidate + 1)
                    obj.population(candidate + 1, :) = trial;
                    obj.energies(candidate + 1) = energy;
                    if energy <= obj.energies(1)
                        obj.promoteLowestEnergy();
                    end
                end
            end
        end
    end
end
