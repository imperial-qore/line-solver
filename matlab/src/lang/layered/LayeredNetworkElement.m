classdef LayeredNetworkElement < Element
    % A generic element of a LayeredNetwork model.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties
        model; % pointer to model
        linConA = []; % Matrix(C,K): admission constraint matrix on this server's layer station
        linConB = []; % Matrix(C,1): admission constraint capacities
        linConRows = {}; % cell of struct('names',cell,'coeffs',vector,'cap',scalar): rows declared by operand name
        lldScaling = []; % vector alpha(n): rate scaling of this server's layer station when it holds n jobs
        lcdScaling = []; % handle beta(n): product-form class-dependent scaling, n counted over this server's operands
        lcdScalingPeak = []; % peak rate scaling per operand, normalizes Util = T*S/peak
        ljdScaling = []; % handle eta(n): non-product-form joint-dependent scaling, n counted over this server's operands
        ljdScalingPeak = []; % peak rate scaling per operand, normalizes Util = T*S/peak
        serverPools = {}; % cell of struct('name','count','rate','compatible'): heterogeneous pools with a compatibility graph over this server's operands
    end

    properties (Hidden, Constant)
        HOST = 0;
        PROCESSOR = 0;
        TASK = 1;
        ENTRY = 2;
        ACTIVITY =3;
        CALL = 4;
    end
    
    
    methods
        %Constructor
        function self = LayeredNetworkElement(name)
            % SELF = LAYEREDNETWORKELEMENT(NAME)
            
            self@Element(name);
        end
        
        function ind = subsindex(self)
            % IND = SUBSINDEX()

            ind = double(self.model.getNodeIndex(self.name))-1 % 0 based
        end

        function self = addConstraint(self, operands, coeffs, cap)
            % self = ADDCONSTRAINT(self, OPERANDS, COEFFS, CAP)
            %
            % Appends one admission constraint row naming its operands, so the
            % meaning does not depend on declaration order:
            %
            %   T2.addConstraint([E2 E3], [1 1], 2);  % n(E2) + n(E3) <= 2
            %   T2.addConstraint(E3, 1, 1);           % n(E3) <= 1
            %
            % OPERANDS are the entries of a Task, or the tasks of a Host, given
            % as handles, names, or a cell mixing the two. COEFFS is a matching
            % vector, or a scalar applied to every operand, and defaults to all
            % ones. Names are resolved against the model in
            % LayeredNetwork.getStruct, where an operand that does not belong to
            % this server is an error rather than a silent mis-mapping.
            if ~isa(self,'Task') && ~isa(self,'Host')
                line_error(mfilename,'Admission constraints can only be set on a Task or a Host, which are the only elements that become server stations in a layer.');
            end
            names = LayeredNetworkElement.operandNames(operands);
            n = length(names);
            if nargin < 3 || isempty(coeffs)
                coeffs = ones(1,n);
            end
            coeffs = double(coeffs(:))';
            if isscalar(coeffs)
                coeffs = coeffs * ones(1,n);
            end
            if length(coeffs) ~= n
                line_error(mfilename,'Admission constraint has %d operands but %d coefficients.', n, length(coeffs));
            end
            if any(~isfinite(coeffs)) || any(coeffs < 0)
                line_error(mfilename,'Admission constraint coefficients must be finite and non-negative.');
            end
            if all(coeffs == 0)
                line_error(mfilename,'Admission constraint has all-zero coefficients, which constrains nothing.');
            end
            if length(unique(names)) ~= n
                line_error(mfilename,'Admission constraint names the same operand more than once; give it a single combined coefficient instead.');
            end
            if nargin < 4 || ~isscalar(cap) || ~isfinite(cap) || cap < 1
                line_error(mfilename,'Admission constraint capacity must be a finite scalar of at least 1.');
            end
            row.names = names;
            row.coeffs = coeffs;
            row.cap = double(cap);
            self.linConRows{end+1,1} = row;
        end

        function self = setConstraint(self, A, b)
            % self = SETCONSTRAINT(self, A, B)
            %
            % Raw form of ADDCONSTRAINT, for programmatic construction. Declares
            % A*n <= B on the station that represents this server in its layer,
            % where n counts the jobs in service or queueing at that station.
            % Columns of A are indexed positionally by the entries of a Task, or
            % by the tasks of a Host, in declaration order, so the mapping shifts
            % if an entry is added later; prefer ADDCONSTRAINT, which names its
            % operands. Only the column count is checked, in
            % LayeredNetwork.getStruct, since entries may be added after this
            % call. Rows from both forms are concatenated.
            if ~isa(self,'Task') && ~isa(self,'Host')
                line_error(mfilename,'Admission constraints can only be set on a Task or a Host, which are the only elements that become server stations in a layer.');
            end
            A = double(A);
            b = double(b(:));
            if isempty(A) || isempty(b)
                line_error(mfilename,'Constraint matrix A and capacity vector b must be non-empty.');
            end
            if ndims(A) > 2 %#ok<ISMAT>
                line_error(mfilename,'Constraint matrix A must be two-dimensional.');
            end
            if size(A,1) ~= length(b)
                line_error(mfilename,'A and b must have matching number of rows.');
            end
            if any(~isfinite(A(:))) || any(A(:) < 0)
                line_error(mfilename,'Constraint matrix A must be finite and non-negative.');
            end
            if any(~isfinite(b)) || any(b < 1)
                line_error(mfilename,'Capacity vector b must be finite and at least 1.');
            end
            if any(all(A == 0, 2))
                line_error(mfilename,'Constraint matrix A has an all-zero row, which constrains nothing.');
            end
            self.linConA = A;
            self.linConB = b;
        end

        function [A, b] = getLinearConstraints(self)
            % [A, B] = GETLINEARCONSTRAINTS(self)

            A = self.linConA;
            b = self.linConB;
        end

        function tf = hasLinearConstraints(self)
            % TF = HASLINEARCONSTRAINTS(self)

            tf = (~isempty(self.linConA) && ~isempty(self.linConB)) || ~isempty(self.linConRows);
        end

        function self = setLoadDependence(self, alpha)
            % self = SETLOADDEPENDENCE(self, ALPHA)
            %
            % ALPHA(n) is the service-rate scaling of the station that
            % represents this server in its layer when that station holds n
            % jobs in total, as in Queue.setLoadDependence. The scaling
            % multiplies the station rate on top of its multiplicity, so a
            % multi-server host applies min(n,m)*ALPHA(n).
            self.assertRateDependent('Load');
            alpha = double(alpha(:))';
            if isempty(alpha) || any(~isfinite(alpha)) || any(alpha <= 0)
                line_error(mfilename,'Load-dependence scalings must be finite and positive.');
            end
            self.lldScaling = alpha;
        end

        function self = setClassDependence(self, beta, peakRatePerOperand)
            % self = SETCLASSDEPENDENCE(self, BETA, PEAKRATEPEROPERAND)
            %
            % BETA(n) is a function handle taking the per-operand population
            % vector of this server: n(j) counts the jobs held on behalf of
            % operand j, which is task j of a Host or entry j of a Task, in the
            % same tasksof/entriesof order as the columns of SETCONSTRAINT. It
            % returns a scalar shared by every operand, or a per-operand vector.
            % PEAKRATEPEROPERAND is REQUIRED (scalar or per-operand vector) and
            % normalizes Util = T*S/peak. Product form holds only where an
            % operand occupies the layer station through a single job class;
            % otherwise SolverLN emits the equivalent joint dependence, which is
            % numerically identical but carries no exactness guarantee.
            self.assertRateDependent('Class');
            if nargin < 3
                peakRatePerOperand = [];
            end
            LayeredNetworkElement.assertDependenceHandle(beta, peakRatePerOperand, 'Class');
            self.lcdScaling = beta;
            self.lcdScalingPeak = double(peakRatePerOperand(:))';
        end

        function self = setJointDependence(self, eta, peakRatePerOperand)
            % self = SETJOINTDEPENDENCE(self, ETA, PEAKRATEPEROPERAND)
            %
            % ETA(n) reads the per-operand population vector of this server
            % arbitrarily (e.g. min(n(1),c)) and is therefore non-product-form:
            % solvers treat it as an approximation. Operand order and the
            % required PEAKRATEPEROPERAND are as in SETCLASSDEPENDENCE.
            self.assertRateDependent('Joint');
            if nargin < 3
                peakRatePerOperand = [];
            end
            LayeredNetworkElement.assertDependenceHandle(eta, peakRatePerOperand, 'Joint');
            if ~isempty(self.serverPools)
                line_error(mfilename,'%s already declares server pools, which are themselves a rate law, so it cannot also take a joint dependence.', self.name);
            end
            self.ljdScaling = eta;
            self.ljdScalingPeak = double(peakRatePerOperand(:))';
        end

        function self = addServerType(self, serverType)
            % self = ADDSERVERTYPE(self, SERVERTYPE)
            %
            % Declares one pool of SERVERTYPE.numOfServers identical servers,
            % each running at SERVERTYPE.rate, eligible only for the operands
            % listed in SERVERTYPE.compatibleClasses:
            %
            %   P1.addServerType(ServerType('Fast', 2, [T2]));
            %   P1.addServerType(ServerType('Shared', 1, [T2 T3]));
            %
            % The operands are the tasks of a Host, or the entries of a Task,
            % given as handles or names. They are resolved against the model in
            % LayeredNetwork.getStruct, where an operand that does not belong to
            % this server is an error rather than a silent mis-mapping, exactly
            % as for ADDCONSTRAINT.
            %
            % SolverLN lowers the whole declaration to the activated-server rate
            % of SN_COMPAT_RATE, carried onto the layer station as a joint
            % dependence, so the pools are an APPROXIMATION in a layer for the
            % same reason SETJOINTDEPENDENCE is.
            self.assertRateDependent('Compatibility');
            if ~isempty(self.ljdScaling)
                line_error(mfilename,'%s already declares a joint dependence, so it cannot also declare server pools, which are a rate law of their own.', self.name);
            end
            if ~isa(serverType,'ServerType')
                line_error(mfilename,'addServerType expects a ServerType object.');
            end
            if serverType.getNumOfServers < 1
                line_error(mfilename,'Server pool ''%s'' must hold at least one server.', serverType.getName);
            end
            compat = serverType.getCompatibleClasses;
            if isempty(compat)
                line_error(mfilename,'Server pool ''%s'' is compatible with no operand, so it can never serve.', serverType.getName);
            end
            % Resolved to names here, inside the class, exactly as ADDCONSTRAINT
            % does: getStruct reads plain strings and cannot reach operandNames.
            compatNames = LayeredNetworkElement.operandNames(compat);
            if numel(unique(compatNames)) ~= numel(compatNames)
                line_error(mfilename,'Server pool ''%s'' names the same operand more than once.', serverType.getName);
            end
            for k = 1:numel(self.serverPools)
                if strcmp(self.serverPools{k}.name, serverType.getName)
                    line_error(mfilename,'Server pool ''%s'' is already declared on %s.', serverType.getName, self.name);
                end
            end
            serverType.setParentQueue(self);
            serverType.setId(numel(self.serverPools));
            pool = struct('name', serverType.getName, ...
                'count', serverType.getNumOfServers, ...
                'rate', serverType.getRate, ...
                'compatible', {compatNames});
            self.serverPools{end+1} = pool;
        end

        function pools = getServerTypes(self)
            % POOLS = GETSERVERTYPES(self) Declared compatibility pools.

            pools = self.serverPools;
        end

        function tf = hasServerPools(self)
            % TF = HASSERVERPOOLS(self)

            tf = ~isempty(self.serverPools);
        end

        function tf = hasRateDependence(self)
            % TF = HASRATEDEPENDENCE(self)

            tf = ~isempty(self.lldScaling) || ~isempty(self.lcdScaling) || ...
                ~isempty(self.ljdScaling) || ~isempty(self.serverPools);
        end

        function [alpha, beta, betaPeak, eta, etaPeak] = getRateDependence(self)
            % [ALPHA, BETA, BETAPEAK, ETA, ETAPEAK] = GETRATEDEPENDENCE(self)

            alpha = self.lldScaling;
            beta = self.lcdScaling;
            betaPeak = self.lcdScalingPeak;
            eta = self.ljdScaling;
            etaPeak = self.ljdScalingPeak;
        end
    end

    methods (Access = private)
        function assertRateDependent(self, what)
            % ASSERTRATEDEPENDENT(self, WHAT) rejects elements that are not layer stations
            if ~isa(self,'Task') && ~isa(self,'Host')
                line_error(mfilename,'%s-dependence can only be set on a Task or a Host, which are the only elements that become server stations in a layer.', what);
            end
            switch SchedStrategy.fromText(self.scheduling)
                case {SchedStrategy.PS, SchedStrategy.FCFS}
                    % the only disciplines whose layer station admits a rate scaling
                otherwise
                    line_error(mfilename,'%s-dependence supported only for processor sharing (PS) and first-come first-serve (FCFS) servers, but %s is scheduled %s.', what, self.name, self.scheduling);
            end
        end
    end

    methods (Static, Access = private)
        function assertDependenceHandle(f, peak, what)
            % ASSERTDEPENDENCEHANDLE(F, PEAK, WHAT) common validation of a dependence declaration
            if ~isa(f,'function_handle')
                line_error(mfilename,'%s dependence must be specified through a function handle.', what);
            end
            if isempty(peak)
                line_error(mfilename,'%s dependence requires an explicit peak rate: pass a scalar (identical peak for every operand) or a per-operand vector.', what);
            end
            if ~isnumeric(peak) || any(~isfinite(peak(:))) || any(peak(:) <= 0)
                line_error(mfilename,'peakRatePerOperand must be a finite positive scalar or per-operand vector.');
            end
        end
    end

    methods (Static, Access = private)
        function names = operandNames(operands)
            % NAMES = OPERANDNAMES(OPERANDS) operand names from handles, names or a cell of both
            if isempty(operands)
                line_error(mfilename,'Admission constraint requires at least one operand.');
            end
            if ischar(operands) || isstring(operands)
                names = {char(operands)};
                return
            end
            if ~iscell(operands)
                operands = num2cell(operands);
            end
            names = cell(1,length(operands));
            for k = 1:length(operands)
                op = operands{k};
                if ischar(op) || isstring(op)
                    names{k} = char(op);
                elseif isa(op,'LayeredNetworkElement')
                    names{k} = op.name;
                else
                    line_error(mfilename,'Admission constraint operand %d is neither a name nor a LayeredNetwork element.', k);
                end
            end
        end
    end
    
end
