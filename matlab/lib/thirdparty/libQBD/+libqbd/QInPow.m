classdef QInPow
    properties
        matrices = {}
        power = 0
        process
    end

    methods
        function obj = QInPow(proc)
            if nargin == 0
                return;
            end

            obj.process = proc;
            obj.power = 1;
            a0 = proc.all_A_0();
            ap = proc.all_A_plus();
            am = proc.all_A_minus();
            if isempty(a0) || isempty(ap) || isempty(am)
                return;
            end

            obj.matrices = cell(1, max([numel(am) + 1, numel(a0), numel(ap)]));
            obj.matrices{1} = {a0{1}, ap{1}};

            im = 1;
            i0 = 2;
            ip = 2;
            for k = 2:numel(obj.matrices)
                obj.matrices{k} = {am{im}, a0{i0}, ap{ip}};
                im = min(im + 1, numel(am));
                i0 = min(i0 + 1, numel(a0));
                ip = min(ip + 1, numel(ap));
            end
        end

        function check(obj)
            if isempty(obj.process) || isempty(obj.process.all_A_0())
                libqbd.raise_error('Infinitesimal generator matrix is empty.');
            end
        end

        function ret = inc_power(obj, step)
            ret = libqbd.QInPow();
            ret.power = obj.power + 1;
            ret.process = obj.process;
            ret.matrices = cell(1, numel(obj.matrices) + 1);

            for k0 = (obj.power + 1):(numel(ret.matrices) - 1)
                true_k0 = min(k0, numel(obj.matrices) - 1);
                ret.matrices{k0 + 1} = {obj.matrices{true_k0 + 1}{1} * obj.process.get_A_minus(k0 - obj.power) * step};
            end

            for k0 = 0:(numel(ret.matrices) - 1)
                true_k0 = min(k0, numel(obj.matrices) - 1);
                blocks = obj.matrices{true_k0 + 1};
                if isempty(ret.matrices{k0 + 1})
                    ret.matrices{k0 + 1} = {};
                end
                for j0 = 0:(numel(blocks) - 1)
                    m = zeros(size(blocks{j0 + 1}));
                    if k0 >= obj.power
                        col = k0 + j0 - obj.power;
                    else
                        col = j0;
                    end

                    left = max(j0 - 1, 0);
                    right = min(j0 + 1, numel(blocks) - 1);
                    if col == 0
                        A = {obj.process.get_A_0(0), obj.process.get_A_minus(1)};
                    else
                        A = {obj.process.get_A_plus(col - 1), obj.process.get_A_0(col), obj.process.get_A_minus(col + 1)};
                    end

                    z_up = max(col - 1, 0);
                    z_left = max(k0 - obj.power, 0);
                    pos = 0;
                    if z_up < z_left && col ~= 0
                        pos = 1;
                    end

                    for i0 = left:right
                        m = m + blocks{i0 + 1} * A{pos + 1} * step;
                        pos = pos + 1;
                    end
                    ret.matrices{k0 + 1}{end + 1} = m;
                end
            end

            for k0 = 0:(numel(ret.matrices) - 1)
                true_k0 = min(k0, numel(obj.matrices) - 1);
                ret.matrices{k0 + 1}{end + 1} = obj.matrices{true_k0 + 1}{end} * obj.process.get_A_plus(k0 + obj.power) * step;
            end
        end

        function obj = mull_by_const(obj, cons)
            for k = 1:numel(obj.matrices)
                for j = 1:numel(obj.matrices{k})
                    obj.matrices{k}{j} = obj.matrices{k}{j} * cons;
                end
            end
        end

        function obj = add_identity_matrix(obj)
            pos0 = 0;
            for k0 = 0:(numel(obj.matrices) - 1)
                m = obj.matrices{k0 + 1}{pos0 + 1};
                obj.matrices{k0 + 1}{pos0 + 1} = m + eye(size(m));
                if k0 < obj.power
                    pos0 = pos0 + 1;
                end
            end
        end

        function obj = plus_assign(obj, right)
            if obj.power >= right.power
                for k0 = 0:(numel(obj.matrices) - 1)
                    rk0 = min(k0, numel(right.matrices) - 1);
                    if k0 < (right.power + 1)
                        first = 0;
                    else
                        first = (numel(obj.matrices{k0 + 1}) - numel(right.matrices{rk0 + 1})) / 2;
                    end
                    for j0 = first:(first + numel(right.matrices{rk0 + 1}) - 1)
                        obj.matrices{k0 + 1}{j0 + 1} = obj.matrices{k0 + 1}{j0 + 1} + right.matrices{rk0 + 1}{j0 - first + 1};
                    end
                end
            else
                tmp = obj;
                obj = right;
                obj = obj.plus_assign(tmp);
            end
        end

        function ret = mull_by_vector(obj, pi)
            ret = cell(1, numel(pi));
            for k = 1:numel(pi)
                ret{k} = zeros(1, numel(pi{k}));
            end

            c = 0;
            for k = numel(pi):(min(numel(obj.matrices), numel(pi) + obj.power) - 1)
                ret{end + 1} = zeros(1, size(obj.matrices{k + 1}{1}, 1));
                c = c + 1;
            end
            for k = 1:(obj.power - c)
                ret{end + 1} = zeros(1, size(obj.matrices{end}{1}, 1));
            end

            matrix_num = 0;
            col_num = 0;
            for k0 = 0:(numel(pi) - 1)
                pi_k = libqbd.as_row_vector(pi{k0 + 1});
                if k0 > obj.power
                    col_num = col_num + 1;
                end
                for j0 = 0:(numel(obj.matrices{matrix_num + 1}) - 1)
                    ret{col_num + j0 + 1} = ret{col_num + j0 + 1} + pi_k * obj.matrices{matrix_num + 1}{j0 + 1};
                end
                if (matrix_num + 1) < numel(obj.matrices)
                    matrix_num = matrix_num + 1;
                end
            end

            while numel(ret) > 1 && max(ret{end}) <= 0.0
                ret(end) = [];
            end
        end
    end
end
