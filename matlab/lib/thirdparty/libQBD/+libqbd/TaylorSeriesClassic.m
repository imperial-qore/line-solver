classdef TaylorSeriesClassic < handle
    properties (Access = private)
        h = 0.0
        B
        is_process_not_binded = true
    end

    methods
        function obj = TaylorSeriesClassic(proc, order, step)
            if nargin > 0
                if nargin < 3
                    step = -1.0;
                end
                obj.bind(proc, order, step);
            end
        end

        function bind(obj, proc, order, step)
            if nargin < 4
                step = -1.0;
            end

            if step < 0
                obj.h = step / proc.get_min_element();
            elseif step > 0
                obj.h = step;
            else
                libqbd.raise_error('step must be not equal 0.');
            end

            obj.B = libqbd.QInPow(proc);
            obj.B = obj.B.mull_by_const(obj.h);
            obj.compute_right_matrix(order);
            obj.is_process_not_binded = false;
        end

        function value = get_step(obj)
            obj.check();
            value = obj.h;
        end

        function dist = get_dist(obj, max_time, pi_0)
            obj.check();
            pi = pi_0;
            dist = {pi};
            time = 0.0;
            while true
                pi = obj.B.mull_by_vector(pi);
                dist{end + 1} = pi;
                time = time + obj.h;
                if time > max_time
                    break;
                end
            end
        end

        function values = get_mean_clients(obj, max_time, pi_0)
            obj.check();
            pi = pi_0;
            values = obj.compute_clients(pi);
            time = 0.0;
            while true
                pi = obj.B.mull_by_vector(pi);
                values(end + 1) = obj.compute_clients(pi); %#ok<AGROW>
                time = time + obj.h;
                if time > max_time
                    break;
                end
            end
        end

        function values = get_mean_queue(obj, queue_size_vector, max_time, pi_0)
            obj.check();
            pi = pi_0;
            values = obj.compute_queue(queue_size_vector, pi);
            time = 0.0;
            while true
                pi = obj.B.mull_by_vector(pi);
                values(end + 1) = obj.compute_queue(queue_size_vector, pi); %#ok<AGROW>
                time = time + obj.h;
                if time > max_time
                    break;
                end
            end
        end
    end

    methods (Access = private)
        function check(obj)
            if obj.is_process_not_binded
                libqbd.raise_error('Not binded to the process.');
            end
            obj.B.check();
        end

        function compute_right_matrix(obj, order)
            order = min(order, libqbd.get_max_factor());
            P = obj.B;
            for k = 1:(order - 1)
                P = P.inc_power(obj.h);
                tmp = P.mull_by_const(libqbd.get_one_div_by_factor(k));
                obj.B = obj.B.plus_assign(tmp);
            end
            obj.B = obj.B.add_identity_matrix();
        end

        function value = compute_clients(~, pi)
            value = 0.0;
            for k = 2:numel(pi)
                value = value + (k - 1) * sum(pi{k});
            end
        end

        function value = compute_queue(~, queue_size_vector, pi)
            value = 0.0;
            qvec = libqbd.as_row_vector(queue_size_vector{1});
            for k = 2:numel(pi)
                if k <= numel(queue_size_vector)
                    qvec = libqbd.as_row_vector(queue_size_vector{k});
                else
                    qvec = qvec + 1.0;
                end
                value = value + sum(libqbd.as_row_vector(pi{k}) .* qvec);
            end
        end
    end
end
