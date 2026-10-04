classdef StationaryDistribution < handle
    properties (Access = private)
        process
        is_binded = false
        rho
        is_rho_computed = false
        R = []
        G = []
        pi_0_c = {}
        mean_cl = []
        is_mean_clients_computed = false
        sum_from_c_to_inf = []
    end

    methods
        function obj = StationaryDistribution(proc)
            if nargin > 0
                obj.bind(proc);
            end
        end

        function bind(obj, proc)
            obj.process = proc;
            obj.is_binded = true;
            obj.is_rho_computed = false;
            obj.is_mean_clients_computed = false;
            obj.R = [];
            obj.G = [];
            obj.pi_0_c = {};
            obj.mean_cl = [];
            obj.sum_from_c_to_inf = [];
        end

        function value = get_rho(obj)
            obj.compute_rho();
            value = obj.rho;
        end

        function matrix = get_R(obj)
            obj.compute_R();
            matrix = obj.R;
        end

        function matrix = get_G(obj)
            obj.compute_G();
            matrix = obj.G;
        end

        function dist = get_dist(obj, max_level)
            obj.compute_pi_0_c();
            dist = obj.pi_0_c(1:min(max_level + 1, numel(obj.pi_0_c)));
            pi = obj.pi_0_c{end};
            for k = numel(dist):max_level
                pi = pi * obj.R;
                dist{end + 1} = pi;
            end
        end

        function dist = get_pi_0_c(obj)
            obj.compute_pi_0_c();
            dist = obj.pi_0_c;
        end

        function value = get_mean_clients(obj)
            if obj.is_mean_clients_computed
                value = obj.mean_cl;
                return;
            end

            obj.compute_rho();
            if obj.rho >= 1
                value = inf;
                return;
            end

            obj.compute_pi_0_c();
            value = 0.0;
            for k = 2:(numel(obj.pi_0_c) - 1)
                value = value + (k - 1) * sum(obj.pi_0_c{k});
            end

            obj.compute_R();
            I = eye(size(obj.R));
            tmp = (I - obj.R) \ I;
            pi = obj.pi_0_c{end};
            value = value + sum(pi * ((obj.R * tmp) + (numel(obj.pi_0_c) - 1) * I) * tmp);
            obj.mean_cl = value;
            obj.is_mean_clients_computed = true;
        end

        function value = get_sum_from_c_to_inf(obj)
            if isempty(obj.sum_from_c_to_inf)
                obj.compute_pi_0_c();
                I = eye(size(obj.R));
                obj.sum_from_c_to_inf = ((I - obj.R).') \ obj.pi_0_c{end}.';
                obj.sum_from_c_to_inf = obj.sum_from_c_to_inf.';
            end
            value = obj.sum_from_c_to_inf;
        end

        function value = get_mean_queue(obj, queue_size_vector)
            obj.compute_pi_0_c();
            if numel(queue_size_vector) ~= numel(obj.pi_0_c)
                libqbd.raise_error('You need to specify c first levels.');
            end

            obj.compute_rho();
            if obj.rho >= 1
                value = inf;
                return;
            end

            value = 0.0;
            for k = 1:(numel(obj.pi_0_c) - 1)
                value = value + sum(libqbd.as_row_vector(obj.pi_0_c{k}) .* libqbd.as_row_vector(queue_size_vector{k}));
            end

            I = eye(size(obj.R));
            tmp = (I - obj.R).';
            pi = obj.pi_0_c{end}.';
            queue_tail = libqbd.as_row_vector(queue_size_vector{end});
            value = value + sum(((tmp \ pi).') .* queue_tail);
            value = value + sum((tmp \ (tmp \ ((obj.pi_0_c{end} * obj.R).'))));
        end
    end

    methods (Access = private)
        function check(obj)
            if ~obj.is_binded
                libqbd.raise_error('Not binded to the process.');
            end
            if isempty(obj.process.all_A_0())
                libqbd.raise_error('Infinitesimal generator matrix is empty.');
            end
        end

        function compute_rho(obj)
            if obj.is_rho_computed
                return;
            end

            obj.check();
            am = obj.process.all_A_minus();
            a0 = obj.process.all_A_0();
            ap = obj.process.all_A_plus();
            A = am{end} + a0{end} + ap{end};
            A(:, 1) = 1.0;
            r = zeros(size(A, 1), 1);
            r(1) = 1.0;
            alpha = (A.' \ r).';
            obj.rho = sum(alpha * ap{end}) / sum(alpha * am{end});
            obj.is_rho_computed = true;
        end

        function compute_G(obj)
            if ~isempty(obj.G)
                return;
            end

            obj.check();
            obj.compute_rho();
            if obj.rho >= 1
                libqbd.raise_error('rho is equal or greater than 1.');
            end

            am = obj.process.all_A_minus();
            a0 = obj.process.all_A_0();
            ap = obj.process.all_A_plus();
            A_m = am{end};
            A_0 = a0{end};
            A_p = ap{end};
            T = -(A_0 \ eye(size(A_0)));
            V_m = T * A_m;
            V_p = T * A_p;
            I = eye(size(T));
            W = (I - V_m * V_p - V_p * V_m) \ I;
            U = I;
            G = V_m;
            tol = realmin(class(G));
            while true
                U = U * V_p;
                V_m = W * V_m * V_m;
                V_p = W * V_p * V_p;
                W = (I - V_m * V_p - V_p * V_m) \ I;
                T_prev = G;
                G = G + U * V_m;
                delta = T_prev - G;
                if max(abs(delta(:))) <= tol
                    break;
                end
            end
            obj.G = G;
        end

        function compute_R(obj)
            if ~isempty(obj.R)
                return;
            end

            obj.check();
            obj.compute_G();
            a0 = obj.process.all_A_0();
            ap = obj.process.all_A_plus();
            A_0 = a0{end};
            A_p = ap{end};
            U = -(A_0 + A_p * obj.G);
            obj.R = A_p / U;
        end

        function compute_pi_0_c(obj)
            obj.check();
            if ~isempty(obj.pi_0_c)
                return;
            end

            obj.compute_rho();
            if obj.rho >= 1
                libqbd.raise_error('rho is equal or greater than 1.');
            end

            c = numel(obj.process.all_A_0());
            if c == 0
                libqbd.raise_error('Generator matrix is empty.');
            elseif c > 1
                c = c - 1;
            end
            c = max(c, 1);

            matrix_len = 0;
            for k = 0:c
                matrix_len = matrix_len + size(obj.process.get_A_0(k), 1);
            end

            B = zeros(matrix_len, matrix_len);
            A00 = obj.process.get_A_0(0);
            Ap0 = obj.process.get_A_plus(0);
            B(1:size(A00, 1), 1:size(A00, 2)) = A00;
            B(1:size(Ap0, 1), (size(A00, 2) + 1):(size(A00, 2) + size(Ap0, 2))) = Ap0;

            row_offset = size(A00, 1) + 1;
            col_offset = 0;
            for k = 1:(c - 1)
                A_minus = obj.process.get_A_minus(k);
                A_0 = obj.process.get_A_0(k);
                A_plus = obj.process.get_A_plus(k);
                B(row_offset:(row_offset + size(A_minus, 1) - 1), (col_offset + 1):(col_offset + size(A_minus, 2))) = A_minus;
                col_offset = col_offset + size(A_minus, 2);
                B(row_offset:(row_offset + size(A_0, 1) - 1), (col_offset + 1):(col_offset + size(A_0, 2))) = A_0;
                B(row_offset:(row_offset + size(A_plus, 1) - 1), (col_offset + size(A_0, 2) + 1):(col_offset + size(A_0, 2) + size(A_plus, 2))) = A_plus;
                row_offset = row_offset + size(A_0, 1);
            end

            A_minus = obj.process.get_A_minus(c);
            B(row_offset:(row_offset + size(A_minus, 1) - 1), (col_offset + 1):(col_offset + size(A_minus, 2))) = A_minus;
            col_offset = col_offset + size(A_minus, 2);
            obj.compute_R();
            A_tail = obj.process.get_A_0(c) + obj.R * obj.process.get_A_minus(c + 1);
            B(row_offset:(row_offset + size(A_tail, 1) - 1), (col_offset + 1):(col_offset + size(A_tail, 2))) = A_tail;

            norm_eq = ones(size(B, 1), 1);
            I = eye(size(obj.R));
            ones_tail = ones(size(obj.R, 1), 1);
            norm_eq(end - size(obj.R, 1) + 1:end) = (I - obj.R) \ ones_tail;
            B(:, 1) = norm_eq;
            right = zeros(size(B, 1), 1);
            right(1) = 1.0;
            dist = (B.' \ right).';
            dist(dist < 0.0) = 0.0;

            left = 1;
            for k = 0:c
                right_idx = left + size(obj.process.get_A_0(k), 1) - 1;
                obj.pi_0_c{end + 1} = dist(left:right_idx);
                left = right_idx + 1;
            end
        end
    end
end
