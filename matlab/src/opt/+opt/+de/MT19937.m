classdef MT19937 < handle
    % MT19937  Bit-exact reimplementation of numpy's legacy MT19937 core,
    % matching the generator underlying numpy.random.RandomState.
    %
    % Reproduces numpy's exact 32-bit output stream so that NumpyRandomState
    % and the differential-evolution optimizer generate identical draws to the
    % native-Python line-opt for a given seed. Seeding follows numpy's
    % _legacy_seeding: a scalar seed that fits in 32 bits uses init_genrand;
    % otherwise the seed words are fed to init_by_array.

    properties
        mt   % 1x624 uint32 state
        mti  % scalar position (numpy get_state pos)
    end

    methods
        function obj = MT19937(seed)
            obj.mt = zeros(1, 624, 'uint32');
            obj.seed(seed);
        end

        function seed(obj, s)
            s = uint64(s);
            if s <= uint64(4294967295)
                obj.initGenrand(uint32(s));
            else
                key = uint32([]);
                while s > 0
                    key(end+1) = uint32(bitand(s, uint64(4294967295))); %#ok<AGROW>
                    s = bitshift(s, -32);
                end
                obj.initByArray(key);
            end
        end

        function initGenrand(obj, s)
            obj.mt(1) = uint32(s);
            for ii = 1:623
                prev = uint64(obj.mt(ii));
                x = bitxor(prev, bitshift(prev, -30));
                val = uint64(1812433253) * x + uint64(ii);
                obj.mt(ii+1) = uint32(bitand(val, uint64(4294967295)));
            end
            obj.mti = 624;
        end

        function initByArray(obj, initKey)
            obj.initGenrand(uint32(19650218));
            nk = numel(initKey);
            % numpy-style 0-based ii,jj; MATLAB mt[ii] == mt(ii+1), mt[ii-1] == mt(ii).
            ii = 1; jj = 0;
            k = max(624, nk);
            while k > 0
                prev = uint64(obj.mt(ii));          % mt[i-1] in numpy == mt(ii) here (ii is 1-based == numpy i)
                mixed = bitxor(prev, bitshift(prev, -30));
                term = uint64(bitand(uint64(mixed) * uint64(1664525), uint64(4294967295)));
                val = bitxor(uint64(obj.mt(ii+1)), term);
                val = mod(uint64(val) + uint64(initKey(jj+1)) + uint64(jj), uint64(4294967296));
                obj.mt(ii+1) = uint32(val);
                ii = ii + 1; jj = jj + 1;
                if ii+1 > 624
                    obj.mt(1) = obj.mt(624);
                    ii = 1;
                end
                if jj >= nk
                    jj = 0;
                end
                k = k - 1;
            end
            for k = 1:623
                prev = uint64(obj.mt(ii));
                mixed = bitxor(prev, bitshift(prev, -30));
                term = uint64(bitand(uint64(mixed) * uint64(1566083941), uint64(4294967295)));
                val = bitxor(uint64(obj.mt(ii+1)), term);
                val = mod(uint64(val) - uint64(ii) + uint64(4294967296), uint64(4294967296));
                obj.mt(ii+1) = uint32(val);
                ii = ii + 1;
                if ii+1 > 624
                    obj.mt(1) = obj.mt(624);
                    ii = 1;
                end
            end
            obj.mt(1) = uint32(2147483648);  % 0x80000000
            obj.mti = 624;
        end

        function y = nextUint32(obj)
            N = 624; M = 397;
            MATRIX_A = uint32(2567483615);   % 0x9908b0df
            UPPER = uint32(2147483648);      % 0x80000000
            LOWER = uint32(2147483647);      % 0x7fffffff
            if obj.mti >= N
                for kk = 1:(N-M)
                    yv = bitor(bitand(obj.mt(kk), UPPER), bitand(obj.mt(kk+1), LOWER));
                    mag = uint32(0);
                    if bitand(yv, uint32(1)) ~= 0
                        mag = MATRIX_A;
                    end
                    obj.mt(kk) = bitxor(bitxor(obj.mt(kk+M), bitshift(yv, -1)), mag);
                end
                for kk = (N-M+1):(N-1)
                    yv = bitor(bitand(obj.mt(kk), UPPER), bitand(obj.mt(kk+1), LOWER));
                    mag = uint32(0);
                    if bitand(yv, uint32(1)) ~= 0
                        mag = MATRIX_A;
                    end
                    obj.mt(kk) = bitxor(bitxor(obj.mt(kk+(M-N)), bitshift(yv, -1)), mag);
                end
                yv = bitor(bitand(obj.mt(N), UPPER), bitand(obj.mt(1), LOWER));
                mag = uint32(0);
                if bitand(yv, uint32(1)) ~= 0
                    mag = MATRIX_A;
                end
                obj.mt(N) = bitxor(bitxor(obj.mt(M), bitshift(yv, -1)), mag);
                obj.mti = 0;
            end

            y = obj.mt(obj.mti + 1);
            obj.mti = obj.mti + 1;
            y = bitxor(y, bitshift(y, -11));
            y = bitxor(y, bitand(bitshift(y, 7), uint32(2636928640)));   % 0x9d2c5680
            y = bitxor(y, bitand(bitshift(y, 15), uint32(4022730752)));  % 0xefc60000
            y = bitxor(y, bitshift(y, -18));
        end

        function key = getStateKey(obj)
            key = obj.mt;
        end
    end
end
