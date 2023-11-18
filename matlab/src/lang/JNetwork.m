classdef JNetwork < Model
    % JLINE extended queueing network model.
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.

    properties (Access=public)
        obj; % java object
    end

    % PUBLIC METHODS
    methods (Access=public)

        %Constructor
        function self = JNetwork(model)
            % SELF = NETWORK(MODELNAME)
            self@Model(model); % model is the model's name
            self.obj = jline.lang.Network(model);
        end

        function sn = getStruct(self, wantInitialState) % get abritrary representation
            if nargin<2
                wantInitialState = false;
            end
            sn = JLINE.from_jline_struct(self.obj, self.obj.getStruct(wantInitialState));
        end

        function R = getNumberOfClasses(self)
            % R = GETNUMBEROFCLASSES()

            R = self.obj.getNumberOfClasses();
        end

        function M = getNumberOfStations(self)
            % M = GETNUMBEROFSTATIONS()

            M = self.obj.getNumberOfStations();
        end


        function I = getNumberOfNodes(self)
            % I = GETNUMBEROFNODES()

            I = self.obj.getNumberOfNodes();
        end

        function bool = hasFork(self)
            % to be changed
            bool = false;
        end

        function ind = getNodeIndex(self, name)
            if ischar(name)
                ind = self.obj.getNodeByName(name);
            else
                ind = self.obj.getNodeIndex(name.obj);
            end
        end

        function self = link(self, P)
            nodes = self.obj.getNodes();
            classes = self.obj.getClasses();
            routing_matrix = jline.lang.RoutingMatrix(self.obj, classes, nodes);
            I = self.obj.getNumberOfNodes;
            R = self.obj.getNumberOfClasses;
            if iscell(P)
                for i = 1:I
                    for j = 1:I
                        for r=1:R
                            for s=1:R
                                if P{r,s}(i,j) > 0
                                    routing_matrix.addConnection(nodes.get(i-1), nodes.get(j-1), classes.get(r-1), classes.get(s-1), P{r,s}(i,j));
                                end
                            end
                        end
                    end
                end
            else % if a single matrix
                for i = 1:I
                    for j = 1:I
                        if P(i,j) > 0
                            routing_matrix.addConnection(nodes.get(i-1), nodes.get(j-1), classes.get(0), P(i,j));
                        end
                    end
                end
            end
            self.obj.link(routing_matrix);
        end

        function P = initRoutingMatrix(self)
            M = self.getNumberOfNodes();
            R = self.getNumberOfClasses();
            P = cellzeros(R,R,M,M);
        end

        function self = setDoChecks(self, bool)
            self.obj.setDoChecks(bool);
        end

    end
end
