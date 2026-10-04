classdef Node < NetworkElement
    % An abstract for a node in a Network model
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties
        model;
        input;
        server;
        output;
        index;
    end

    methods(Hidden)
        %Constructor
        function self = Node(name)
            % SELF = NODE(NAME)

            self@NetworkElement(char(name));
            self.index = NaN;
        end

        function self = setModel(self, model)
            % SELF = SETMODEL(MODEL)
            %
            % Add a pointer to the model object

            self.model = model;
        end

        function self = link(self, nodeTo)
            % SELF = LINK(NODETO)
            %
            %

            self.model.addLink(self,nodeTo);
        end

        function self = reset(self)
            % SELF = RESET()
            %
            % Reset internal data structures when the network model is
            % reset

        end
    end

    methods

        function sections = getSections(self)
            % SECTIONS = GETSECTIONS()

            sections = {self.input, self.server, self.output};
        end

        function remaining = remainingClassIndexes(self, jobclass)
            % REMAINING = REMAININGCLASSINDEXES(JOBCLASS)
            %
            % Indexes of the classes that survive the removal of JOBCLASS.
            % The lookup is by name so that a class object belonging to
            % another copy of the model still resolves, as in
            % @MNetwork/removeClass.m.

            K = length(self.model.classes);
            r = self.model.getClassByName(jobclass.name).index;
            remaining = setdiff(1:K, r);
        end

        function self = removeJobClass(self, jobclass)
            % SELF = REMOVEJOBCLASS(JOBCLASS)
            %
            % Remove all per-class configuration referencing JOBCLASS from
            % this node. The base implementation drops the class's routing
            % (output) strategy; subclasses extend it to drop service,
            % capacity, arrival and class-switching configuration. Called by
            % @MNetwork/removeClass.m, mirroring Node.removeJobClass in the
            % JAR and Node.remove_job_class in python.

            remaining = self.remainingClassIndexes(jobclass);
            if ~isempty(self.output) && isprop(self.output, 'outputStrategy') ...
                    && numel(self.output.outputStrategy) == numel(remaining) + 1
                self.output.outputStrategy = self.output.outputStrategy(remaining);
            end
        end

        function setProbRouting(self, class, destination, probability)
            % SETPROBROUTING(CLASS, DESTINATION, PROBABILITY)

            setRouting(self, class, RoutingStrategy.PROB, destination, probability);
        end

        function setRouting(self, class, strategy, par1, par2)
            % SETROUTING(CLASS, STRATEGY, PARAM)
            % SETROUTING(CLASS, STRATEGY, DESTINATION, PROBABILITY)

            %global GlobalConstants.CoarseTol

            if self.model.isJavaNative()
                jline_classes = self.model.obj.getClasses();
                switch strategy
                    case RoutingStrategy.RAND
                        self.obj.setRouting(jline_classes.get(class.index-1),jline.lang.constant.RoutingStrategy.RAND);
                    case RoutingStrategy.RROBIN
                        self.obj.setRouting(jline_classes.get(class.index-1),jline.lang.constant.RoutingStrategy.RROBIN);
                    case RoutingStrategy.WRROBIN
                        node_target =  self.model.obj.getNodeByName(par1.getName());
                        weight = par2;
                        self.obj.setRouting(jline_classes.get(class.index-1),jline.lang.constant.RoutingStrategy.WRROBIN, node_target, weight);
                    case RoutingStrategy.DISABLED
                        self.obj.setRouting(jline_classes.get(class.index-1),jline.lang.constant.RoutingStrategy.DISABLED);
                    case RoutingStrategy.PROB
                        line_error(mfilename, 'Use setProbRouting to assign routing probabilities for a JNetwork node.');
                end
                return
            end
            
            if isa(self,'Cache')
                switch strategy
                    case {RoutingStrategy.SQ, RoutingStrategy.WRROBIN, RoutingStrategy.RROBIN}
                        line_error(mfilename,'State-dependent routing not supported with caches. Add instead a Router node after the cache.');
                end
            end
            switch strategy
                case RoutingStrategy.SQ
                    % SETROUTING(CLASS, RoutingStrategy.SQ, D) - SQ(d): sample D
                    % destinations uniformly with replacement and route to the
                    % shortest of them. Dispatcher memory is not supported.
                    if nargin < 4
                        par1 = 2; % d: sampled destinations
                    end
                    param_d = par1;
                    if ~isscalar(param_d) || param_d < 1 || abs(param_d-round(param_d)) > GlobalConstants.CoarseTol
                        line_error(mfilename,'SQ parameter d must be a positive integer.');
                    end
                    if nargin >= 5 && ~isempty(par2) && any(double(par2) ~= 0)
                        line_error(mfilename,'SQ with dispatcher memory is not supported. Only SQ(d) is available.');
                    end
                    self.output.outputStrategy{1, class.index}{2} = RoutingStrategy.toText(strategy);
                    self.output.outputStrategy{1, class.index}{3}{1} = round(param_d);
                case RoutingStrategy.WRROBIN                    
                    destination = par1;
                    weight = par2;
                    if abs(weight-round(weight)) < GlobalConstants.CoarseTol
                        self.output.outputStrategy{1, class.index}{2} = RoutingStrategy.toText(strategy);
                        if length(self.output.outputStrategy{1, class.index})<3
                            self.output.outputStrategy{1, class.index}{3}{1} = {destination, weight};
                        else
                            pos = Node.findRoutingEntry(self.output.outputStrategy{1, class.index}{3}, destination);
                            if pos > 0
                                self.output.outputStrategy{1, class.index}{3}{pos} = {destination, weight};
                            else
                                self.output.outputStrategy{1, class.index}{3}{end+1} = {destination, weight};
                            end
                        end
                    else
                        line_error(mfilename,'Weighted round robin weights must be integers.')
                    end
                otherwise
                    switch nargin
                        case 3 % no destination specified
                            self.output.outputStrategy{1, class.index}{2} = RoutingStrategy.toText(strategy);
                        case 5
                            destination = par1;
                            probability = par2;
                            self.output.outputStrategy{1, class.index}{2} = RoutingStrategy.toText(strategy);
                            if length(self.output.outputStrategy{1, class.index})<3
                                self.output.outputStrategy{1, class.index}{3}{1} = {destination, probability};
                            else
                                % re-declaring (class,destination) overwrites, never appends
                                pos = Node.findRoutingEntry(self.output.outputStrategy{1, class.index}{3}, destination);
                                if pos > 0
                                    self.output.outputStrategy{1, class.index}{3}{pos} = {destination, probability};
                                else
                                    self.output.outputStrategy{1, class.index}{3}{end+1} = {destination, probability};
                                end
                            end
                    end
            end
        end

        function setStateDepRouting(self, class, departure, branches, level, C, d)
            % SETSTATEDEPROUTING(CLASS, DEPARTURE, BRANCHES, LEVEL, C, D)
            %
            % Declares this node to be the entry center e of a subnetwork
            % Q(V,V) served by the product-form state-dependent routing of
            % Krzesinski (1987), "Multiclass Queueing Networks with
            % State-Dependent Routing", Performance Evaluation 7:125-143.
            %
            % DEPARTURE is the departure center d of Q(V,V), which may be this
            % node itself in a central server model. BRANCHES is a cell array
            % following the paper's own indexing: BRANCHES{1} must be empty
            % because branch index 1 denotes the complement M-V, and
            % BRANCHES{b} for b >= 2 lists the nodes of branch b with its entry
            % center first and its departure center last. A single-center
            % branch is written {node}. LEVEL(b) is the index t of the
            % subnetwork with B_b in V_t - V_{t+1}, and LEVEL(1) is ignored. C
            % is the 1xT vector of coefficients C_t and D the TxB matrix of
            % coefficients d_tb, of which entry (t,b) is read for
            % 1 <= t <= LEVEL(b) and 2 <= b <= B.
            %
            % Negative C_t and positive d_tb make the routing prefer the least
            % congested branches and impose the population bounds
            % m_b <= d_{t,b}/(-C_t) and v_t <= D_tt/(-C_t).
            %
            % Example, the central server of Section 2.5 with two peripheral
            % centers, C = (-1,-1), d_12 = 1, d_13 = 2, d_23 = 3:
            %   d = zeros(2,3); d(1,2) = 1; d(1,3) = 2; d(2,3) = 3;
            %   node1.setStateDepRouting(class, node1, {[], {node2}, {node3}}, ...
            %                            [0 1 2], [-1 -1], d);

            if ~iscell(branches) || isempty(branches) || ~isempty(branches{1})
                line_error(mfilename,'BRANCHES must be a cell array whose first entry is empty: branch index 1 denotes the complement M-V.');
            end
            B = numel(branches);
            T = numel(C);
            if numel(level) ~= B
                line_error(mfilename,'LEVEL must have one entry per branch index, including the unused index 1.');
            end
            if size(d,1) < T || size(d,2) < B
                line_error(mfilename,sprintf('D must be at least %dx%d.',T,B));
            end
            sdr = struct();
            sdr.departure = departure;
            sdr.branch = cell(1,B);
            sdr.entryOf = cell(1,B);
            sdr.departureOf = cell(1,B);
            for b = 2:B
                bnodes = branches{b};
                if ~iscell(bnodes)
                    bnodes = {bnodes};
                end
                if isempty(bnodes)
                    line_error(mfilename,sprintf('Branch %d is empty.',b));
                end
                sdr.branch{b} = bnodes;
                sdr.entryOf{b} = bnodes{1};
                sdr.departureOf{b} = bnodes{end};
            end
            sdr.level = level(:)';
            sdr.C = C(:)';
            sdr.d = d;
            self.output.outputStrategy{1, class.index}{2} = RoutingStrategy.toText(RoutingStrategy.SDR);
            self.output.outputStrategy{1, class.index}{3} = {sdr};
        end

        function bool = hasClassSwitching(self)
            % BOOL = HASCLASSSWITCHING()

            bool = isa(self.server,'ClassSwitcher');
        end

        function bool = isStateful(self)
            % BOOL = ISSTATEFUL()

            bool = isa(self,'StatefulNode');
        end

        function bool = isStation(self)
            % BOOL = ISSTATION()

            bool = isa(self,'Station');
        end
    end

    methods(Access = protected)
        % Override copyElement method:
        function clone = copyElement(self)
            % CLONE = COPYELEMENT()

            % Make a shallow copy of all properties
            clone = copyElement@Copyable(self);
            % Make a deep copy of each object
            clone.input = self.input.copy;
            clone.server = self.server.copy;
            clone.output = self.output.copy;
        end
    end

    methods (Access = public)
        function ind = subsindex(self)
            % IND = SUBSINDEX()
            if isa(self.model,'Network')
                % Handle the new delegation pattern
                node_idx = self.model.getNodeIndex(self.name);
                if isempty(node_idx) || isnan(node_idx) || node_idx <= 0
                    error('Node:subsindex', 'Invalid node index for %s: %g', self.name, node_idx);
                end
                ind = double(node_idx)-1; % 0 based for MATLAB indexing
            elseif isa(self.model,'MNetwork') 
                node_idx = self.model.getNodeIndex(self.name);
                if isempty(node_idx) || isnan(node_idx) || node_idx <= 0
                    error('Node:subsindex', 'Invalid node index for %s: %g', self.name, node_idx);
                end
                ind = double(node_idx)-1; % 0 based
            elseif isa(self.model,'JNetwork') 
                ind = self.model.obj.getNodeIndex(self.obj);
            else
                error('Node:subsindex', 'Unsupported model type: %s', class(self.model));
            end
        end

        function V = horzcat(self, varargin)
            % V = HORZCAT(VARARGIN)

            V = zeros(1, length(varargin) + 1);
            try
                self_idx = self.subsindex();
                if numel(self_idx) ~= 1
                    error('Node:horzcat', 'subsindex returned non-scalar value for %s: %s', self.name, mat2str(self_idx));
                end
                V(1) = 1+ self_idx;
                
                for v=1:length(varargin)
                    if isa(varargin{v}, 'Node')
                        node_idx = varargin{v}.subsindex();
                        if numel(node_idx) ~= 1
                            error('Node:horzcat', 'subsindex returned non-scalar value for %s: %s', varargin{v}.name, mat2str(node_idx));
                        end
                        V(1+v) = 1+ node_idx;
                    else
                        error('Node:horzcat', 'Element %d is not a Node object', v);
                    end
                end
            catch e
                error('Node:horzcat', 'Error in horizontal concatenation: %s', e.message);
            end
        end

        function V = vertcat(self, varargin)
            % V = VERTCAT(VARARGIN)

            V = zeros(length(varargin) + 1, 1);
            try
                V(1) = 1+ self.subsindex;
                for v=1:length(varargin)
                    if isa(varargin{v}, 'Node')
                        V(1+v) = 1+varargin{v}.subsindex;
                    else
                        error('Node:vertcat', 'Element %d is not a Node object', v);
                    end
                end
            catch e
                error('Node:vertcat', 'Error in vertical concatenation: %s', e.message);
            end
        end

        function summary(self)
            % SUMMARY()

            line_printf('\nNode: <strong>%s</strong>',self.getName);
            %self.input.summary;
            %            self.server.summary;
            %            self.output.summary;
        end
    end

    methods (Static)
        function pos = findRoutingEntry(entries, destination)
            % POS = FINDROUTINGENTRY(ENTRIES, DESTINATION)
            % Index of the {destination, value} pair naming DESTINATION, 0 if absent.
            pos = 0;
            if isempty(entries) || ~isa(destination,'Node')
                return
            end
            for e = 1:length(entries)
                entry = entries{e};
                if iscell(entry) && ~isempty(entry) && isa(entry{1},'Node') && strcmp(entry{1}.name, destination.name)
                    pos = e;
                    return
                end
            end
        end
    end
end
