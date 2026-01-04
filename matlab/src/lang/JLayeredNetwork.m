classdef JLayeredNetwork < Model
    % JLINE extended layered queueing network model.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Access=public)
        obj; % java object
    end

    % PUBLIC METHODS
    methods (Access=public)

        %Constructor
        function self = JLayeredNetwork(name)
            % SELF = JLAYEREDNETWORK(MODELNAME)
            self@Model(name); % model is the model's name
            self.obj = jline.lang.layered.LayeredNetwork(name);
        end

        function lqn = getStruct(self, wantInitialState) % get arbitrary representation
            if nargin<2
                wantInitialState = false;
            end
            lqn = JLINE.from_jline_layered_struct(self.obj, self.obj.getStruct(wantInitialState));
        end

        function E = getNumberOfLayers(self)
            % E = GETNUMBEROFLAYERS()
            E = self.obj.getNumberOfLayers();
        end

        function E = getNumberOfModels(self)
            % E = GETNUMBEROFMODELS()
            E = self.obj.getNumberOfModels();
        end

        function H = getNumberOfHosts(self)
            % H = GETNUMBEROFHOSTS()
            H = self.obj.getNumberOfHosts();
        end

        function T = getNumberOfTasks(self)
            % T = GETNUMBEROFTASKS()
            T = self.obj.getNumberOfTasks();
        end

        function E = getNumberOfEntries(self)
            % E = GETNUMBEROFENTRIES()
            E = self.obj.getNumberOfEntries();
        end

        function A = getNumberOfActivities(self)
            % A = GETNUMBEROFACTIVITIES()
            A = self.obj.getNumberOfActivities();
        end

        function ind = getNodeIndex(self, name)
            if ischar(name)
                ind = self.obj.getNodeByName(name);
            else
                ind = self.obj.getNodeIndex(name.obj);
            end
        end

        function node = getNodeByName(self, name)
            % NODE = GETNODEBYNAME(NAME)
            
            javaNode = self.obj.getNodeByName(name);
            if isempty(javaNode)
                node = NaN;
            else
                % Create a wrapper MATLAB object for the Java node
                % We need to determine the node type and create appropriate wrapper
                nodeClassName = char(javaNode.getClass().getSimpleName());
                switch nodeClassName
                    case 'Host'
                        node = Host(self, char(javaNode.getName()));
                        node.obj = javaNode;
                    case 'Processor'
                        node = Processor(self, char(javaNode.getName()));
                        node.obj = javaNode;
                    case 'Task'
                        node = Task(self, char(javaNode.getName()));
                        node.obj = javaNode;
                    case 'CacheTask'
                        node = CacheTask(self, char(javaNode.getName()));
                        node.obj = javaNode;
                    case 'FunctionTask'
                        node = FunctionTask(self, char(javaNode.getName()));
                        node.obj = javaNode;
                    case 'Entry'
                        node = Entry(self, char(javaNode.getName()));
                        node.obj = javaNode;
                    case 'ItemEntry'
                        node = ItemEntry(self, char(javaNode.getName()));
                        node.obj = javaNode;
                    case 'Activity'
                        node = Activity(self, char(javaNode.getName()));
                        node.obj = javaNode;
                    otherwise
                        % Fallback - create generic LayeredNetworkElement wrapper
                        node = LayeredNetworkElement(char(javaNode.getName()));
                        node.obj = javaNode;
                        node.model = self;
                end
            end
        end

        function self = addHost(self, host)
            % SELF = ADDHOST(HOST)
            self.obj.addHost(host.obj);
        end

        function self = addTask(self, task)
            % SELF = ADDTASK(TASK)
            self.obj.addTask(task.obj);
        end

        function self = addEntry(self, entry)
            % SELF = ADDENTRY(ENTRY)
            self.obj.addEntry(entry.obj);
        end

        function self = addActivity(self, activity)
            % SELF = ADDACTIVITY(ACTIVITY)
            self.obj.addActivity(activity.obj);
        end

        function reset(self)
            self.obj.reset();
        end

        function sanitize(self)
            % SANITIZE()
            % Validates the LayeredNetwork configuration.
            % Ensures that if entries are defined, activities are also defined to serve those entries.
            self.obj.sanitize();
        end

        function writeXML(self, filename)
            % WRITEXML(FILENAME)
            self.obj.writeXML(filename);
        end

    end

    methods (Static)

        function model = parseXML(filename)
            % MODEL = PARSEXML(FILENAME)
            javaModel = jline.lang.layered.LayeredNetwork.parseXML(filename);
            model = JLayeredNetwork(char(javaModel.getName()));
            model.obj = javaModel;
        end

    end
end