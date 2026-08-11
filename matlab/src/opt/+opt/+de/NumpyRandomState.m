classdef NumpyRandomState < handle
    % NumpyRandomState  Bit-exact port of the subset of numpy.random.RandomState
    % consumed by scipy's differential_evolution. Backed by opt.de.MT19937.
    %
    % Implements random_sample, uniform, the default-dtype randint (32-bit
    % masked rejection), and shuffle/permutation (Fisher-Yates driven by
    % random_interval), each reproducing numpy's exact draw sequence and word
    % consumption, so the optimizer follows the identical trajectory to the
    % native-Python line-opt for a given integer seed.

    properties
        gen  % opt.de.MT19937 handle
    end

    methods
        function obj = NumpyRandomState(seed)
            obj.gen = opt.de.MT19937(seed);
        end

        function d = randomSample(obj)
            a = double(bitshift(obj.gen.nextUint32(), -5));   % 27 bits
            b = double(bitshift(obj.gen.nextUint32(), -6));   % 26 bits
            d = (a * 67108864 + b) / 9007199254740992;
        end

        function out = randomSampleN(obj, n)
            out = zeros(1, n);
            for i = 1:n
                out(i) = obj.randomSample();
            end
        end

        function d = uniformScalar(obj, low, high)
            d = low + (high - low) * obj.randomSample();
        end

        function out = uniformN(obj, low, high, n)
            out = zeros(1, n);
            for i = 1:n
                out(i) = low + (high - low) * obj.randomSample();
            end
        end

        function m = fillMask(~, v)
            v = uint64(v);
            v = bitor(v, bitshift(v, -1));
            v = bitor(v, bitshift(v, -2));
            v = bitor(v, bitshift(v, -4));
            v = bitor(v, bitshift(v, -8));
            v = bitor(v, bitshift(v, -16));
            m = v;
        end

        function v = randint(obj, low, high)
            rng_ = uint64(high) - uint64(1) - uint64(low);
            if rng_ == 0
                v = low;
                return;
            end
            mask = obj.fillMask(rng_);
            while true
                w = uint64(obj.gen.nextUint32());
                val = bitand(w, mask);
                if val <= rng_
                    v = double(low) + double(val);
                    return;
                end
            end
        end

        function v = randomInterval(obj, maxv)
            maxv = uint64(maxv);
            if maxv == 0
                v = 0;
                return;
            end
            mask = obj.fillMask(maxv);
            while true
                w = uint64(obj.gen.nextUint32());
                val = bitand(w, mask);
                if val <= maxv
                    v = double(val);
                    return;
                end
            end
        end

        function arr = shuffle(obj, arr)
            n = numel(arr);
            for p = n:-1:2
                jj = obj.randomInterval(uint64(p - 1));   % 0-based j in [0, p-1]
                tmp = arr(p);
                arr(p) = arr(jj + 1);
                arr(jj + 1) = tmp;
            end
        end

        function arr = permutation(obj, n)
            arr = 0:(n - 1);
            arr = obj.shuffle(arr);
        end
    end
end
