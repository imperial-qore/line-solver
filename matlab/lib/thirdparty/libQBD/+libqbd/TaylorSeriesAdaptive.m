classdef TaylorSeriesAdaptive < handle
    properties (Access = private)
        is_binded = false
        max_saved_data = 0
        max_degree = 0
        proc
        dist_in_ref_points = {}
        norms_in_ref_points = {}
        derivs_in_ref_points = {}
        ref_points = []
        clients_in_ref_points = []
        min_elem = 0.0
        error = 0.0
    end

    methods
        function obj = TaylorSeriesAdaptive(proc, pi0, error_value, max_time, approx_type, strategy)
            if nargin > 0
                if nargin < 5
                    approx_type = libqbd.APPROX_TAYLOR_UNLIM();
                end
                if nargin < 6
                    strategy = libqbd.STRATEGY_FAST();
                end
                obj.bind(proc, pi0, error_value, max_time, approx_type, strategy);
            end
        end

        function bind(obj, proc, pi0, error_value, max_time, approx_type, strategy)
            if nargin < 6
                approx_type = libqbd.APPROX_TAYLOR_UNLIM();
            end
            if nargin < 7
                strategy = libqbd.STRATEGY_FAST();
            end
            if obj.is_binded
                libqbd.raise_error('Already binded.');
            end

            obj.min_elem = -proc.get_min_element();
            obj.proc = proc;
            obj.dist_in_ref_points = {pi0};
            obj.norms_in_ref_points = {{}};
            obj.derivs_in_ref_points = {{}};
            obj.ref_points = 0.0;
            obj.error = error_value;
            obj.max_degree = libqbd.parse_approx_type(approx_type);
            obj.max_saved_data = min(obj.max_degree, double(strategy));
            while obj.ref_points(end) < max_time
                obj.next_point(libqbd.get_max_factor());
            end
            obj.is_binded = true;
        end

        function values = get_reference_times(obj)
            obj.check();
            values = obj.ref_points;
        end

        function values = get_reference_dists(obj)
            obj.check();
            values = obj.dist_in_ref_points;
        end

        function [dist, errors] = get_dist(obj, times)
            obj.check();
            errors = [];
            if isempty(times)
                dist = {};
                return;
            end

            while obj.ref_points(end) < times(end)
                obj.next_point(libqbd.get_max_factor());
            end

            num = obj.find_max_le_t(times(1), 1);
            norm_times = [];
            dist = {};
            errors = [];
            for t = times
                reg = (t - obj.ref_points(num)) * obj.min_elem;
                if reg > 1.0
                    [dist, errors] = obj.calc_intermed_points(dist, errors, norm_times, num);
                    norm_times = [];
                    num = obj.find_max_le_t(t, num);
                    reg = (t - obj.ref_points(num)) * obj.min_elem;
                end
                norm_times(end + 1) = reg; %#ok<AGROW>
            end

            if ~isempty(norm_times)
                [dist, errors] = obj.calc_intermed_points(dist, errors, norm_times, num);
            end
        end

        function values = get_reference_mean_clients(obj)
            obj.check();
            if numel(obj.clients_in_ref_points) < numel(obj.ref_points)
                for k = (numel(obj.clients_in_ref_points) + 1):numel(obj.ref_points)
                    obj.clients_in_ref_points(end + 1) = libqbd.function_of_dist(obj.dist_in_ref_points{k}, @client_weights); %#ok<AGROW>
                end
            end
            values = obj.clients_in_ref_points;

            function weights = client_weights(level, length)
                weights = ones(1, length) * level;
            end
        end
    end

    methods (Access = private)
        function check(obj)
            if ~obj.is_binded
                libqbd.raise_error('Not binded to the process.');
            end
            if isempty(obj.proc.all_A_0())
                libqbd.raise_error('Infinitesimal generator matrix is empty.');
            end
        end

        function next_point(obj, n)
            deriv = obj.dist_in_ref_points{end};
            res = deriv;
            norms = {};
            derivs = {};
            k = 0;
            two_delta_inv = 0.5;
            two_delta_in_n = two_delta_inv;
            min_elem_inv = 1.0 / obj.min_elem;
            er = inf;
            while er > obj.error && k < n
                deriv = obj.proc.mull_by_row_vector(deriv, min_elem_inv);
                res = libqbd.vec_fma(res, deriv, libqbd.get_one_div_by_factor(k));
                norm_value = libqbd.l1norm_cell(deriv);
                er = norm_value * gammainc(2.0, k + 2.0, 'lower') * exp(2.0) * two_delta_in_n;
                if numel(norms) < obj.max_saved_data
                    norms{end + 1} = norm_value; %#ok<AGROW>
                    derivs{end + 1} = deriv; %#ok<AGROW>
                end
                two_delta_in_n = two_delta_in_n * two_delta_inv;
                k = k + 1;
            end
            obj.dist_in_ref_points{end + 1} = res;
            obj.ref_points(end + 1) = obj.ref_points(end) + 1.0 / obj.min_elem;
            obj.norms_in_ref_points{end + 1} = norms;
            obj.derivs_in_ref_points{end + 1} = derivs;
        end

        function [dist, errors] = calc_intermed_points(obj, dist, errors, deltas, num)
            for d = deltas
                if d <= eps
                    dist{end + 1} = obj.dist_in_ref_points{num}; %#ok<AGROW>
                    errors(end + 1) = 0.0; %#ok<AGROW>
                    continue;
                end

                res = obj.dist_in_ref_points{num};
                deriv_comp = {};
                k = 0;
                delta_m = d;
                delta_inv = 0.5;
                delta_minus_n = delta_inv;
                er = inf;
                while k < obj.max_degree && er > obj.error
                    if k < numel(obj.norms_in_ref_points{num})
                        deriv = obj.derivs_in_ref_points{num}{k + 1};
                        norm_value = obj.norms_in_ref_points{num}{k + 1};
                    else
                        if isempty(deriv_comp)
                            if ~isempty(obj.derivs_in_ref_points{num})
                                deriv_comp = obj.derivs_in_ref_points{num}{end};
                            else
                                deriv_comp = obj.dist_in_ref_points{num};
                            end
                        end
                        deriv_comp = obj.proc.mull_by_row_vector(deriv_comp, 1.0 / obj.min_elem);
                        deriv = deriv_comp;
                        norm_value = libqbd.l1norm_cell(deriv_comp);
                    end

                    res = libqbd.vec_fma(res, deriv, delta_m * libqbd.get_one_div_by_factor(k));
                    delta_m = delta_m * d;
                    er = norm_value * gammainc(2.0 * d, k + 2.0, 'lower') * exp(2.0 * d) * delta_minus_n;
                    delta_minus_n = delta_minus_n * delta_inv;
                    k = k + 1;
                end
                dist{end + 1} = res; %#ok<AGROW>
                errors(end + 1) = er; %#ok<AGROW>
            end
        end

        function num = find_max_le_t(obj, t, num)
            idx = find(obj.ref_points(num:end) <= t, 1, 'last');
            if isempty(idx)
                num = 1;
            else
                num = num + idx - 1;
            end
        end
    end
end
