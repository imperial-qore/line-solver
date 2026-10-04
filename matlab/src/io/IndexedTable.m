classdef IndexedTable < handle
    % INDEXEDTABLE Generic enhanced table wrapper with flexible object-based filtering
    %
    % This class wraps any MATLAB table to enable filtering using Station, Node,
    % JobClass, and/or Chain objects, while maintaining full backward compatibility
    % with standard table operations.
    %
    % Supports all NetworkSolver result tables:
    % - AvgTable (Station + JobClass)
    % - AvgNodeTable (Node + JobClass)
    % - AvgChainTable (Station + Chain)
    % - AvgNodeChainTable (Node + Chain)
    % - AvgSysTable (Chain only)
    % - Any custom table with Station/Node/JobClass/Chain columns
    %
    % Usage:
    %   % Create from any MATLAB table
    %   t_native = solver.avgTable();
    %   t = IndexedTable(t_native);
    %
    %   % Standard table operations (works as normal)
    %   rows = t.data(1:5, :)
    %   value = t.data.QLen(1)
    %
    %   % Object-based filtering (four equivalent syntaxes)
    %   rows = t(station_obj)                      % Direct indexing
    %   rows = t.filterBy(station_obj)             % filterBy method
    %   rows = t.get(station_obj)                  % get alias (shorthand)
    %   rows = t.tget(station_obj)                 % tget alias
    %
    %   % Filter by class object
    %   rows = t(class_obj)
    %   rows = t.filterBy(class_obj)
    %   rows = t.get(class_obj)
    %   rows = t.tget(class_obj)
    %
    %   % Filter by both station/node and class/chain (any order)
    %   rows = t(station_obj, class_obj)
    %   rows = t.filterBy(station_obj, class_obj)
    %   rows = t.get(station_obj, class_obj)
    %   rows = t.tget(station_obj, class_obj)
    %
    % Example (AvgTable):
    %   solver = SolverJMT(model);
    %   solver.runAnalyzer();
    %   avgTable = IndexedTable(solver.avgTable());
    %   metrics = avgTable(queue, jobclass);
    %
    % Example (AvgNodeTable):
    %   nodeTable = IndexedTable(solver.avgNodeTable());
    %   metrics = nodeTable(sink, jobclass);
    %
    % Example (AvgChainTable):
    %   chainTable = IndexedTable(solver.avgChainTable());
    %   metrics = chainTable(queue, chain);
    %
    % Example (AvgSysTable):
    %   sysTable = IndexedTable(solver.avgSysTable());
    %   metrics = sysTable(chain);

    properties
        data  % Underlying MATLAB table
    end

    methods
        function obj = IndexedTable(table_data)
            % INDEXEDTABLE(table_data) Create wrapper around any MATLAB table
            %
            % Args:
            %   table_data: A MATLAB table object with object-filterable columns

            if ~istable(table_data)
                error('IndexedTable requires a MATLAB table as input');
            end
            obj.data = table_data;
        end

        function varargout = subsref(obj, s)
            % SUBSREF Support both standard indexing and object-based filtering
            %
            % Intelligently detects object types and filters appropriate columns

            if isempty(s)
                varargout = {obj};
                return;
            end

            switch s(1).type
                case '()'
                    % Check if this is object-based indexing
                    indices = s(1).subs;

                    % Check if first argument is a Station/Node/JobClass/Chain object
                    isObjArg = isa(indices{1}, 'Station') || isa(indices{1}, 'Node') || ...
                               isa(indices{1}, 'JobClass') || isa(indices{1}, 'Chain');

                    if isObjArg
                        % Object-based indexing - use filterBy
                        if length(indices) == 1
                            result = obj.filterBy(indices{1});
                        elseif length(indices) == 2
                            result = obj.filterBy(indices{1}, indices{2});
                        else
                            error('Too many arguments for object-based indexing');
                        end
                    else
                        % Standard numeric/logical indexing - pass to underlying table
                        varargout{1} = subsref(obj.data, s);
                        return;
                    end

                    % Handle chained indexing
                    if length(s) > 1
                        varargout{1} = subsref(result, s(2:end));
                        return;
                    end

                    varargout{1} = result;

                case '.'
                    % Handle method calls and property access
                    if ~isempty(s(1).subs)
                        methodName = s(1).subs;
                        % Check if this is a call to our alias methods
                        if strcmp(methodName, 'get') && length(s) > 1 && strcmp(s(2).type, '()')
                            result = obj.get(s(2).subs{:});
                            if length(s) > 2
                                varargout{1} = subsref(result, s(3:end));
                                return;
                            end
                            varargout{1} = result;
                            return;
                        elseif strcmp(methodName, 'tget') && length(s) > 1 && strcmp(s(2).type, '()')
                            result = obj.tget(s(2).subs{:});
                            if length(s) > 2
                                varargout{1} = subsref(result, s(3:end));
                                return;
                            end
                            varargout{1} = result;
                            return;
                        elseif strcmp(methodName, 'filterBy') && length(s) > 1 && strcmp(s(2).type, '()')
                            result = obj.filterBy(s(2).subs{:});
                            if length(s) > 2
                                varargout{1} = subsref(result, s(3:end));
                                return;
                            end
                            varargout{1} = result;
                            return;
                        end
                        % Check if this is a table column name - delegate to underlying data
                        if ismember(methodName, obj.data.Properties.VariableNames)
                            varargout{1} = subsref(obj.data, s);
                            return;
                        end
                    end
                    % For other properties/methods, delegate to builtin subsref
                    try
                        varargout{1} = builtin('subsref', obj, s);
                    catch ME
                        if strcmp(ME.identifier, 'MATLAB:maxlhs')
                            builtin('subsref', obj, s);
                            varargout = {};
                        else
                            rethrow(ME);
                        end
                    end

                case '{}'
                    % Cell indexing - delegate to underlying table
                    varargout{1} = subsref(obj.data, s);

                otherwise
                    error('Invalid indexing type');
            end
        end

        function result = filterBy(obj, varargin)
            % FILTERBY Filter table by Station/Node and/or JobClass/Chain objects
            %
            % Intelligently filters based on table structure and object types
            %
            % Usage:
            %   rows = t.filterBy(station_obj)
            %   rows = t.filterBy(class_obj)
            %   rows = t.filterBy(station_obj, class_obj)
            %
            % Args:
            %   varargin: 1 or 2 arguments - Station/Node/JobClass/Chain objects
            %
            % Returns:
            %   result: Filtered MATLAB table

            if nargin == 2
                % Single argument
                arg1 = varargin{1};
                if isa(arg1, 'Station') || isa(arg1, 'Node')
                    result = obj.filterByFirstDimension(arg1);
                elseif isa(arg1, 'JobClass') || isa(arg1, 'Chain')
                    result = obj.filterBySecondDimension(arg1);
                else
                    error('Invalid argument type. Expected Station, Node, JobClass, or Chain.');
                end

            elseif nargin == 3
                % Two arguments
                arg1 = varargin{1};
                arg2 = varargin{2};

                result = obj.data;

                % Determine order: could be (station, class) or (class, station)
                if (isa(arg1, 'Station') || isa(arg1, 'Node')) && ...
                   (isa(arg2, 'JobClass') || isa(arg2, 'Chain'))
                    % First arg is station/node, second is class/chain
                    result = obj.filterByFirstDimensionInternal(result, arg1);
                    result = obj.filterBySecondDimensionInternal(result, arg2);

                elseif (isa(arg1, 'JobClass') || isa(arg1, 'Chain')) && ...
                       (isa(arg2, 'Station') || isa(arg2, 'Node'))
                    % First arg is class/chain, second is station/node
                    result = obj.filterByFirstDimensionInternal(result, arg2);
                    result = obj.filterBySecondDimensionInternal(result, arg1);

                else
                    error('Expected one Station/Node and one JobClass/Chain object.');
                end

            else
                error('filterBy expects 1 or 2 arguments');
            end
        end

        function result = tget(obj, varargin)
            % TGET Alias for filterBy - backward compatibility
            result = obj.filterBy(varargin{:});
        end

        function result = get(obj, varargin)
            % GET Alias for filterBy - convenient shorthand
            result = obj.filterBy(varargin{:});
        end

        function n = numel(obj, varargin)
            % NUMEL Support numel() function
            % When called with extra arguments (during indexing like obj(a,b)),
            % return 1 to tell MATLAB we expect a single output from subsref.
            % When called directly as numel(obj), return actual element count.
            if nargin == 1
                n = numel(obj.data);
            else
                n = 1;
            end
        end

        function n = numArgumentsFromSubscript(obj, s, indexingContext)
            % NUMARGUMENTSFROMSUBSCRIPT Specify number of outputs for indexing
            % This tells MATLAB to always expect exactly 1 output from subsref
            % operations, which is appropriate for table data access.
            n = 1;
        end

        function n = size(obj, varargin)
            % SIZE Support size() function
            n = size(obj.data, varargin{:});
        end

        function n = height(obj)
            % HEIGHT Get number of rows
            n = height(obj.data);
        end

        function n = width(obj)
            % WIDTH Get number of columns
            n = width(obj.data);
        end

        function disp(obj)
            % DISP Display the table
            disp(obj.data);
        end

        function summary(obj)
            % SUMMARY Display table summary
            summary(obj.data);
        end

        function t = head(obj, n)
            % HEAD Get first n rows
            if nargin < 2
                n = 5;
            end
            t = head(obj.data, n);
        end

        function t = tail(obj, n)
            % TAIL Get last n rows
            if nargin < 2
                n = 5;
            end
            t = tail(obj.data, n);
        end

        function t = getTable(obj)
            % GETTABLE Return the native MATLAB table
            t = obj.data;
        end

        function result = varfun(func, obj, varargin)
            % VARFUN Apply function to table variables
            % Delegates to the underlying MATLAB table
            result = varfun(func, obj.data, varargin{:});
        end

        function result = table2array(obj)
            % TABLE2ARRAY Convert IndexedTable to array
            % Delegates to the underlying MATLAB table
            result = table2array(obj.data);
        end
    end

    methods (Access = private)
        function result = filterByFirstDimension(obj, first_dim_obj)
            % Filter table by first dimension object (Station or Node)
            result = obj.filterByFirstDimensionInternal(obj.data, first_dim_obj);
        end

        function result = filterByFirstDimensionInternal(~, tbl, first_dim_obj)
            % Internal method for filtering by first dimension
            colNames = tbl.Properties.VariableNames;

            % Try Station first, then Node
            if ismember('Station', colNames)
                col = tbl.Station;
                if iscategorical(col)
                    result = tbl(col == first_dim_obj.name, :);
                else
                    result = tbl(strcmp(col, first_dim_obj.name), :);
                end
            elseif ismember('Node', colNames)
                col = tbl.Node;
                if iscategorical(col)
                    result = tbl(col == first_dim_obj.name, :);
                else
                    result = tbl(strcmp(col, first_dim_obj.name), :);
                end
            else
                result = tbl;
            end
        end

        function result = filterBySecondDimension(obj, second_dim_obj)
            % Filter table by second dimension object (JobClass or Chain)
            result = obj.filterBySecondDimensionInternal(obj.data, second_dim_obj);
        end

        function result = filterBySecondDimensionInternal(~, tbl, second_dim_obj)
            % Internal method for filtering by second dimension
            colNames = tbl.Properties.VariableNames;

            % Try JobClass first, then Chain
            if ismember('JobClass', colNames)
                col = tbl.JobClass;
                if iscategorical(col)
                    result = tbl(col == second_dim_obj.name, :);
                else
                    result = tbl(strcmp(col, second_dim_obj.name), :);
                end
            elseif ismember('Chain', colNames)
                col = tbl.Chain;
                if iscategorical(col)
                    result = tbl(col == second_dim_obj.name, :);
                else
                    result = tbl(strcmp(col, second_dim_obj.name), :);
                end
            else
                result = tbl;
            end
        end
    end
end
