classdef QBD < handle
    properties (Access = private)
        A_plus = {};
        A_0 = {};
        A_minus = {};
    end

    methods
        function values = all_A_plus(obj)
            values = obj.A_plus;
        end

        function values = all_A_0(obj)
            values = obj.A_0;
        end

        function values = all_A_minus(obj)
            values = obj.A_minus;
        end

        function matrix = get_A_minus(obj, level)
            if level == 0
                libqbd.raise_error('Matrix A_minus for zero level is undefined.');
            end
            if isempty(obj.A_minus)
                libqbd.raise_error('No levels specified.');
            end

            idx = min(level, numel(obj.A_minus));
            matrix = obj.A_minus{idx};
        end

        function matrix = get_A_0(obj, level)
            if isempty(obj.A_0)
                libqbd.raise_error('No levels specified.');
            end

            idx = min(level + 1, numel(obj.A_0));
            matrix = obj.A_0{idx};
        end

        function matrix = get_A_plus(obj, level)
            if isempty(obj.A_plus)
                libqbd.raise_error('No levels specified.');
            end

            idx = min(level + 1, numel(obj.A_plus));
            matrix = obj.A_plus{idx};
        end

        function add_zero_level(obj, varargin)
            if ~(isempty(obj.A_0) && isempty(obj.A_plus))
                libqbd.raise_error('Level zero already exists.');
            end

            if nargin == 2
                A_plus = varargin{1};
                obj.A_plus{1} = A_plus;
                obj.A_0{1} = diag(-sum(A_plus, 2));
                return;
            end

            A_0 = varargin{1};
            A_plus = varargin{2};
            if size(A_0, 1) ~= size(A_plus, 1)
                libqbd.raise_error('Different number of rows in matrices of the same level.');
            end
            if size(A_0, 1) ~= size(A_0, 2)
                libqbd.raise_error('Matrix A(0) is not square.');
            end

            obj.A_0{1} = A_0;
            obj.A_plus{1} = A_plus;
        end

        function add_level(obj, varargin)
            obj.check_filled_levels();

            if nargin == 3
                A_minus = varargin{1};
                A_plus = varargin{2};
                if size(A_minus, 1) ~= size(A_plus, 1)
                    libqbd.raise_error('Different number of rows in matrices of the same level.');
                end
                if size(A_minus, 2) ~= size(obj.A_0{end}, 2)
                    libqbd.raise_error('The number of columns of the matrix A(-) is not equal to the number of columns of the matrix A(0) of the previous level.');
                end

                obj.A_minus{end + 1} = A_minus;
                obj.A_0{end + 1} = diag(-(sum(A_minus, 2) + sum(A_plus, 2)));
                obj.A_plus{end + 1} = A_plus;
                return;
            end

            A_minus = varargin{1};
            A_0 = varargin{2};
            A_plus = varargin{3};
            if size(A_minus, 1) ~= size(A_0, 1) || size(A_minus, 1) ~= size(A_plus, 1)
                libqbd.raise_error('Different number of rows in matrices of the same level.');
            end
            if size(A_minus, 2) ~= size(obj.A_0{end}, 2)
                libqbd.raise_error('The number of columns of the matrix A(-) is not equal to the number of columns of the matrix A(0) of the previous level.');
            end
            if size(A_0, 2) ~= size(obj.A_plus{end}, 2)
                libqbd.raise_error('The number of columns of the matrix A(0) is not equal to the number of columns of the matrix A(+) of the previous level.');
            end
            if size(A_0, 1) ~= size(A_0, 2)
                libqbd.raise_error('Matrix A(0) is not square.');
            end

            obj.A_minus{end + 1} = A_minus;
            obj.A_0{end + 1} = A_0;
            obj.A_plus{end + 1} = A_plus;
        end

        function add_final_level(obj, varargin)
            obj.check_filled_levels();

            if nargin == 2
                A_minus = varargin{1};
                if size(A_minus, 1) ~= size(obj.A_plus{end}, 1)
                    libqbd.raise_error('Different number of rows in matrices of the same level.');
                end
                if size(A_minus, 2) ~= size(obj.A_0{end}, 2)
                    libqbd.raise_error('The number of columns of the matrix A(-) is not equal to the number of columns of the matrix A(0) of the previous level.');
                end
                if size(A_minus, 1) ~= size(A_minus, 2)
                    libqbd.raise_error('Matrix A(-) is not square.');
                end
                if size(obj.A_plus{end}, 1) ~= size(obj.A_plus{end}, 2)
                    libqbd.raise_error('Matrix A(+) is not square.');
                end

                obj.A_minus{end + 1} = A_minus;
                obj.A_0{end + 1} = diag(-(sum(A_minus, 2) + sum(obj.A_plus{end}, 2)));
                obj.A_plus{end + 1} = obj.A_plus{end};
                return;
            end

            A_minus = varargin{1};
            A_0 = varargin{2};
            if size(A_minus, 1) ~= size(A_0, 1) || size(A_minus, 1) ~= size(obj.A_plus{end}, 1)
                libqbd.raise_error('Different number of rows in matrices of the same level.');
            end
            if size(A_minus, 2) ~= size(obj.A_0{end}, 2)
                libqbd.raise_error('The number of columns of the matrix A(-) is not equal to the number of columns of the matrix A(0) of the previous level.');
            end
            if size(A_0, 2) ~= size(obj.A_plus{end}, 2)
                libqbd.raise_error('The number of columns of the matrix A(0) is not equal to the number of columns of the matrix A(+) of the previous level.');
            end
            if size(A_minus, 1) ~= size(A_minus, 2)
                libqbd.raise_error('Matrix A(-) is not square.');
            end
            if size(A_0, 1) ~= size(A_0, 2)
                libqbd.raise_error('Matrix A(0) is not square.');
            end
            if size(obj.A_plus{end}, 1) ~= size(obj.A_plus{end}, 2)
                libqbd.raise_error('Matrix A(+) is not square.');
            end

            obj.A_minus{end + 1} = A_minus;
            obj.A_0{end + 1} = A_0;
            obj.A_plus{end + 1} = obj.A_plus{end};
        end

        function fix_diagonal(obj)
            if ~(isempty(obj.A_0) || isempty(obj.A_plus))
                obj.A_0{1} = obj.adjust_diagonal(obj.A_0{1}, sum(obj.A_plus{1}, 2) + sum(obj.A_0{1}, 2));
            end

            n = min([numel(obj.A_0), numel(obj.A_plus), numel(obj.A_minus) + 1]);
            for k = 2:n
                delta = sum(obj.A_minus{k - 1}, 2) + sum(obj.A_plus{k}, 2) + sum(obj.A_0{k}, 2);
                obj.A_0{k} = obj.adjust_diagonal(obj.A_0{k}, delta);
            end
        end

        function value = get_min_element(obj)
            value = 0.0;
            for k = 1:numel(obj.A_0)
                diag_values = diag(obj.A_0{k});
                value = min(value, min(diag_values));
            end
        end

        function res = mull_by_row_vector(obj, vec, cons)
            if nargin < 3
                cons = 1.0;
            end

            n = numel(vec);
            if n == 0
                libqbd.raise_error('An empty vector was passed.');
            end

            vec = cellfun(@libqbd.as_row_vector, vec, 'UniformOutput', false);
            res = {};

            if n >= 3
                res{1} = (vec{1} * obj.get_A_0(0) + vec{2} * obj.get_A_minus(1)) * cons;
                for p = 3:n
                    res{p - 1} = (vec{p - 2} * obj.get_A_plus(p - 3) + vec{p - 1} * obj.get_A_0(p - 2) + vec{p} * obj.get_A_minus(p - 1)) * cons;
                end
                res{n} = (vec{n - 1} * obj.get_A_plus(n - 2) + vec{n} * obj.get_A_0(n - 1)) * cons;
                tmp = (vec{n} * obj.get_A_plus(n - 1)) * cons;
                if norm(tmp, 1) > 0.0
                    res{n + 1} = tmp;
                end
            elseif n == 2
                res{1} = (vec{1} * obj.get_A_0(0) + vec{2} * obj.get_A_minus(1)) * cons;
                res{2} = (vec{1} * obj.get_A_plus(0) + vec{2} * obj.get_A_0(1)) * cons;
                res{3} = (vec{2} * obj.get_A_plus(1)) * cons;
            else
                res{1} = (vec{1} * obj.get_A_0(0)) * cons;
                res{2} = (vec{1} * obj.get_A_plus(0)) * cons;
            end
        end
    end

    methods (Access = private)
        function check_filled_levels(obj)
            if numel(obj.A_0) ~= numel(obj.A_plus) || numel(obj.A_plus) ~= (numel(obj.A_minus) + 1)
                libqbd.raise_error('Unfilled levels found.');
            end
        end

        function matrix = adjust_diagonal(~, matrix, delta)
            idx = 1:(size(matrix, 1) + 1):numel(matrix);
            matrix(idx) = matrix(idx) - reshape(delta, 1, []);
        end
    end
end
