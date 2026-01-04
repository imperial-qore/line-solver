classdef NetworkXMLIO < handle
    % NetworkXMLIO - MATLAB wrapper for XML import/export of Network models
    %
    % This class provides MATLAB interface to save and load Network models 
    % to/from XML format using the Java NetworkXMLIO implementation.
    % Supports class switching through ClassSwitch nodes and cross-class routing.
    %
    % Usage:
    %   % Export a network to XML
    %   NetworkXMLIO.exportToXML(model, 'mymodel.xml');
    %
    %   % Import a network from XML  
    %   model = NetworkXMLIO.importFromXML('mymodel.xml');
    %
    %   % Create class switching example
    %   NetworkXMLIO.demonstrateClassSwitching();
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    methods (Static)
        
        function exportToXML(network, filename)
            % EXPORTTOXML(NETWORK, FILENAME)
            %
            % Exports a Network model to XML file
            %
            % Parameters:
            %   network - Network model to export
            %   filename - Output XML file path
            %
            % Example:
            %   model = Network('MyModel');
            %   source = Source(model, 'Source');
            %   queue = Queue(model, 'Queue', SchedStrategy.FCFS);
            %   sink = Sink(model, 'Sink');
            %   NetworkXMLIO.exportToXML(model, 'mymodel.xml');
            
            try
                if isa(network, 'Network')
                    if network.isJavaNative()
                        % Direct Java Network export
                        jline.io.NetworkXMLIO.exportToXML(network.implementation.obj, filename);
                    else
                        % Convert MNetwork to JNetwork for export
                        tempJNetwork = NetworkXMLIO.convertToJNetwork(network);
                        jline.io.NetworkXMLIO.exportToXML(tempJNetwork.obj, filename);
                    end
                else
                    error('NetworkXMLIO:InvalidNetwork', 'Input must be a Network object');
                end
                
                fprintf('Network successfully exported to: %s\n', filename);
                
            catch e
                error('NetworkXMLIO:ExportFailed', 'Failed to export network: %s', e.message);
            end
        end
        
        function network = importFromXML(filename)
            % IMPORTFROMXML(FILENAME)
            %
            % Imports a Network model from XML file
            %
            % Parameters:
            %   filename - Input XML file path
            %
            % Returns:
            %   network - Loaded Network model
            %
            % Example:
            %   model = NetworkXMLIO.importFromXML('mymodel.xml');
            
            try
                % Import using Java implementation
                javaNetwork = jline.io.NetworkXMLIO.importFromXML(filename);
                
                % Wrap in MATLAB Network with Java implementation
                network = Network(char(javaNetwork.getName()), 'java');
                network.implementation.obj = javaNetwork;
                
                fprintf('Network successfully imported from: %s\n', filename);
                
            catch e
                error('NetworkXMLIO:ImportFailed', 'Failed to import network: %s', e.message);
            end
        end
        
        function validateXML(filename)
            % VALIDATEXML(FILENAME)
            %
            % Validates XML file structure without full import
            %
            % Parameters:
            %   filename - XML file path to validate
            
            try
                % Basic validation by attempting to parse the XML
                docFactory = javax.xml.parsers.DocumentBuilderFactory.newInstance();
                docBuilder = docFactory.newDocumentBuilder();
                doc = docBuilder.parse(filename);
                doc.getDocumentElement().normalize();
                
                root = doc.getDocumentElement();
                if ~strcmp(char(root.getNodeName()), 'network')
                    error('NetworkXMLIO:InvalidXML', 'Root element must be <network>');
                end
                
                fprintf('XML file is valid: %s\n', filename);
                
            catch e
                error('NetworkXMLIO:ValidationFailed', 'XML validation failed: %s', e.message);
            end
        end
        
        function info = getXMLInfo(filename)
            % GETXMLINFO(FILENAME)
            %
            % Gets basic information about XML file without full import
            %
            % Parameters:
            %   filename - XML file path
            %
            % Returns:
            %   info - Struct with network information
            
            try
                docFactory = javax.xml.parsers.DocumentBuilderFactory.newInstance();
                docBuilder = docFactory.newDocumentBuilder();
                doc = docBuilder.parse(filename);
                doc.getDocumentElement().normalize();
                
                root = doc.getDocumentElement();
                
                info.name = char(root.getAttribute('name'));
                info.version = char(root.getAttribute('version'));
                
                % Count nodes
                nodeList = doc.getElementsByTagName('node');
                info.nodeCount = nodeList.getLength();
                
                % Count classes
                classList = doc.getElementsByTagName('class');
                info.classCount = classList.getLength();
                
                % Count connections
                linkList = doc.getElementsByTagName('link');
                info.connectionCount = linkList.getLength();
                
                % Get node types
                info.nodeTypes = {};
                for i = 0:nodeList.getLength()-1
                    nodeElement = nodeList.item(i);
                    nodeType = char(nodeElement.getAttribute('type'));
                    if ~any(strcmp(info.nodeTypes, nodeType))
                        info.nodeTypes{end+1} = nodeType;
                    end
                end
                
                % Count routing matrices (including cross-class routing)
                matrixList = doc.getElementsByTagName('matrix');
                info.routingMatrixCount = matrixList.getLength();
                
                % Check for class switching features
                info.hasClassSwitching = false;
                info.classSwitchNodes = {};
                
                % Check for ClassSwitch nodes
                for i = 0:nodeList.getLength()-1
                    nodeElement = nodeList.item(i);
                    nodeType = char(nodeElement.getAttribute('type'));
                    if strcmp(nodeType, 'ClassSwitch')
                        info.hasClassSwitching = true;
                        nodeName = char(nodeElement.getAttribute('name'));
                        info.classSwitchNodes{end+1} = nodeName;
                    end
                end
                
                % Check for cross-class routing matrices
                info.hasCrossClassRouting = false;
                for i = 0:matrixList.getLength()-1
                    matrixElement = matrixList.item(i);
                    fromClass = char(matrixElement.getAttribute('fromClass'));
                    toClass = char(matrixElement.getAttribute('toClass'));
                    if ~isempty(fromClass) && ~isempty(toClass) && ~strcmp(fromClass, toClass)
                        info.hasCrossClassRouting = true;
                        info.hasClassSwitching = true;
                        break;
                    end
                end
                
            catch e
                error('NetworkXMLIO:InfoFailed', 'Failed to get XML info: %s', e.message);
            end
        end
        
        function demonstrateUsage()
            % DEMONSTRATEUSAGE()
            %
            % Demonstrates basic usage of NetworkXMLIO with a simple example
            
            fprintf('=== NetworkXMLIO Usage Demonstration ===\n');
            
            % Create a simple network
            model = Network('DemoModel');
            
            source = Source(model, 'Source');
            queue = Queue(model, 'Queue', SchedStrategy.FCFS);
            sink = Sink(model, 'Sink');
            
            jobClass = OpenClass(model, 'Class1');
            source.setArrival(jobClass, Exp(1));
            queue.setService(jobClass, Exp(2));
            
            model.link(Network.serialRouting(source, queue, sink));
            
            fprintf('Created demo network with %d nodes\n', model.getNumberOfNodes());
            
            % Export to XML
            filename = 'demo_network.xml';
            NetworkXMLIO.exportToXML(model, filename);
            
            % Show XML info
            info = NetworkXMLIO.getXMLInfo(filename);
            fprintf('XML Info:\n');
            fprintf('  Name: %s\n', info.name);
            fprintf('  Nodes: %d\n', info.nodeCount);
            fprintf('  Classes: %d\n', info.classCount);
            fprintf('  Node Types: %s\n', strjoin(info.nodeTypes, ', '));
            
            % Import back
            importedModel = NetworkXMLIO.importFromXML(filename);
            fprintf('Imported network with %d nodes\n', importedModel.getNumberOfNodes());
            
            fprintf('=== Demonstration Complete ===\n');
        end
        
        function demonstrateClassSwitching()
            % DEMONSTRATECLASSSWITCHING()
            %
            % Demonstrates class switching functionality with NetworkXMLIO
            % Creates a model similar to example_classSwitching_1.m
            
            fprintf('=== Class Switching Demonstration ===\n');
            
            try
                % Create network with class switching (similar to example_classSwitching_1)
                model = Network('ClassSwitchingDemo');
                
                % Create nodes
                source = Source(model, 'Source');
                queue = Queue(model, 'Queue', SchedStrategy.FCFS);
                sink = Sink(model, 'Sink');
                classSwitch = ClassSwitch(model, 'ClassSwitch');
                
                % Create job classes
                class1 = OpenClass(model, 'Class1');
                class2 = OpenClass(model, 'Class2');
                
                % Set arrivals and services
                source.setArrival(class1, Exp.fitMean(10));
                source.setArrival(class2, Exp.fitMean(2));
                queue.setService(class1, Exp.fitMean(1));
                queue.setService(class2, Exp.fitMean(1));
                
                % Set up class switching matrix
                csmatrix = classSwitch.initClassSwitchMatrix();
                csmatrix(class1, class1) = 0.3;
                csmatrix(class1, class2) = 0.7;
                csmatrix(class2, class1) = 1.0;
                classSwitch.setClassSwitchingMatrix(csmatrix);
                
                % Set up routing
                P = model.initRoutingMatrix();
                P{1,1}(source, classSwitch) = 1;
                P{1,1}(classSwitch, queue) = 1;
                P{1,1}(queue, sink) = 1;
                P{2,2}(source, classSwitch) = 1;
                P{2,2}(classSwitch, queue) = 1;
                P{2,2}(queue, sink) = 1;
                model.link(P);
                
                fprintf('Created class switching network with:\n');
                fprintf('  - %d nodes\n', model.getNumberOfNodes());
                fprintf('  - %d classes\n', model.getNumberOfClasses());
                fprintf('  - ClassSwitch node with switching matrix\n');
                
                % Export to XML
                filename = 'class_switching_demo.xml';
                NetworkXMLIO.exportToXML(model, filename);
                
                % Show detailed XML info
                info = NetworkXMLIO.getXMLInfo(filename);
                fprintf('\nXML Analysis:\n');
                fprintf('  Model: %s\n', info.name);
                fprintf('  Nodes: %d (%s)\n', info.nodeCount, strjoin(info.nodeTypes, ', '));
                fprintf('  Classes: %d\n', info.classCount);
                fprintf('  Routing matrices: %d\n', info.routingMatrixCount);
                fprintf('  Has class switching: %s\n', NetworkXMLIO.yesNo(info.hasClassSwitching));
                
                if info.hasClassSwitching
                    if ~isempty(info.classSwitchNodes)
                        fprintf('  ClassSwitch nodes: %s\n', strjoin(info.classSwitchNodes, ', '));
                    end
                    if info.hasCrossClassRouting
                        fprintf('  Has cross-class routing: Yes\n');
                    end
                end
                
                % Import back and validate
                importedModel = NetworkXMLIO.importFromXML(filename);
                fprintf('\nImported network with %d nodes\n', importedModel.getNumberOfNodes());
                
                % Compare structure
                if importedModel.getNumberOfNodes() == model.getNumberOfNodes() && ...
                   importedModel.getNumberOfClasses() == model.getNumberOfClasses()
                    fprintf('✓ Import successful - structure matches\n');
                else
                    fprintf('⚠ Import completed but structure differs\n');
                end
                
            catch e
                fprintf('Class switching demonstration failed: %s\n', e.message);
                fprintf('Note: Class switching requires full MNetwork to JNetwork conversion\n');
            end
            
            fprintf('\n=== Class Switching Demonstration Complete ===\n');
        end
        
        function demonstrateEmbeddedClassSwitching()
            % DEMONSTRATEEMBEDDEDCLASSSWITCHING()
            %
            % Demonstrates embedded class switching in routing matrices
            % Creates a model similar to example_classSwitching_2.m
            
            fprintf('=== Embedded Class Switching Demonstration ===\n');
            
            try
                % Create network with embedded class switching
                model = Network('EmbeddedClassSwitchingDemo');
                
                % Create nodes
                source = Source(model, 'Source');
                queue0 = Queue(model, 'Queue0', SchedStrategy.FCFS);
                queue1 = Queue(model, 'Queue1', SchedStrategy.FCFS);
                queue2 = Queue(model, 'Queue2', SchedStrategy.FCFS);
                sink = Sink(model, 'Sink');
                
                % Create job classes
                class1 = OpenClass(model, 'Class1');
                class2 = OpenClass(model, 'Class2');
                class3 = OpenClass(model, 'Class3');
                
                % Set arrivals and services
                source.setArrival(class1, Exp.fitMean(1));
                queue0.setService(class1, Exp.fitMean(10));
                queue1.setService(class2, Exp.fitMean(20));
                queue2.setService(class3, Exp.fitMean(30));
                
                % Set up cross-class routing matrix (embedded class switching)
                P = model.initRoutingMatrix();
                P{1,1}(source, queue0) = 1;
                P{1,1}(queue0, queue0) = 0.2;  % Stay in same class
                P{1,2}(queue0, queue1) = 0.3;  % Switch to class 2
                P{1,3}(queue0, queue2) = 0.5;  % Switch to class 3
                P{2,2}(queue1, sink) = 1;
                P{3,3}(queue2, sink) = 1;
                model.link(P);
                
                fprintf('Created embedded class switching network with:\n');
                fprintf('  - %d nodes\n', model.getNumberOfNodes());
                fprintf('  - %d classes\n', model.getNumberOfClasses());
                fprintf('  - Cross-class routing probabilities\n');
                
                % Export to XML
                filename = 'embedded_class_switching_demo.xml';
                NetworkXMLIO.exportToXML(model, filename);
                
                % Show detailed XML info
                info = NetworkXMLIO.getXMLInfo(filename);
                fprintf('\nXML Analysis:\n');
                fprintf('  Has cross-class routing: %s\n', NetworkXMLIO.yesNo(info.hasCrossClassRouting));
                fprintf('  Total routing matrices: %d\n', info.routingMatrixCount);
                
            catch e
                fprintf('Embedded class switching demonstration failed: %s\n', e.message);
            end
            
            fprintf('\n=== Embedded Class Switching Demonstration Complete ===\n');
        end
        
        function analyzeClassSwitchingXML(filename)
            % ANALYZECLASSSWITCHINGXML(FILENAME)
            %
            % Analyzes an XML file for class switching features
            %
            % Parameters:
            %   filename - XML file to analyze
            
            if ~exist(filename, 'file')
                error('NetworkXMLIO:FileNotFound', 'File not found: %s', filename);
            end
            
            try
                info = NetworkXMLIO.getXMLInfo(filename);
                
                fprintf('=== Class Switching Analysis: %s ===\n', filename);
                fprintf('Model: %s (version %s)\n', info.name, info.version);
                fprintf('Nodes: %d\n', info.nodeCount);
                fprintf('Classes: %d\n', info.classCount);
                fprintf('Routing matrices: %d\n', info.routingMatrixCount);
                fprintf('\n');
                
                fprintf('Class Switching Features:\n');
                fprintf('  Has class switching: %s\n', NetworkXMLIO.yesNo(info.hasClassSwitching));
                
                if info.hasClassSwitching
                    if ~isempty(info.classSwitchNodes)
                        fprintf('  ClassSwitch nodes: %s\n', strjoin(info.classSwitchNodes, ', '));
                    end
                    fprintf('  Has cross-class routing: %s\n', NetworkXMLIO.yesNo(info.hasCrossClassRouting));
                end
                
                fprintf('  Node types: %s\n', strjoin(info.nodeTypes, ', '));
                
            catch e
                error('NetworkXMLIO:AnalysisFailed', 'Analysis failed: %s', e.message);
            end
        end
        
    end
    
    methods (Static, Access = private)
        
        function jnetwork = convertToJNetwork(mnetwork)
            % CONVERTTOJNETWORK(MNETWORK)
            %
            % Converts MNetwork to JNetwork for XML export
            % Enhanced to handle class switching features
            
            jnetwork = JNetwork(mnetwork.getName());
            
            try
                % Convert nodes
                NetworkXMLIO.convertNodes(mnetwork, jnetwork);
                
                % Convert classes
                NetworkXMLIO.convertClasses(mnetwork, jnetwork);
                
                % Convert routing matrices (including class switching)
                NetworkXMLIO.convertRoutingMatrices(mnetwork, jnetwork);
                
                % Convert service distributions
                NetworkXMLIO.convertServiceDistributions(mnetwork, jnetwork);
                
            catch e
                warning('NetworkXMLIO:ConversionPartial', ...
                       'Partial conversion completed. Some features may not be preserved: %s', e.message);
            end
            
            warning('NetworkXMLIO:ConversionLimited', ...
                   'MNetwork to JNetwork conversion is simplified. Some features may not be preserved.');
        end
        
        function convertNodes(mnetwork, jnetwork)
            % CONVERTNODES(MNETWORK, JNETWORK)
            %
            % Converts nodes from MNetwork to JNetwork
            
            nodes = mnetwork.getNodes();
            for i = 1:length(nodes)
                node = nodes{i};
                nodeName = node.getName();
                
                % Create corresponding Java node based on type
                if isa(node, 'Source')
                    jnode = jline.lang.nodes.Source(jnetwork, nodeName);
                elseif isa(node, 'Queue')
                    jnode = jline.lang.nodes.Queue(jnetwork, nodeName, node.getSchedStrategy());
                    if ~isinf(node.getNumberOfServers())
                        jnode.setNumberOfServers(node.getNumberOfServers());
                    end
                elseif isa(node, 'Sink')
                    jnode = jline.lang.nodes.Sink(jnetwork, nodeName);
                elseif isa(node, 'Router')
                    jnode = jline.lang.nodes.Router(jnetwork, nodeName);
                elseif isa(node, 'ClassSwitch')
                    jnode = jline.lang.nodes.ClassSwitch(jnetwork, nodeName);
                    % Note: Class switching matrix conversion would go here
                elseif isa(node, 'Delay')
                    jnode = jline.lang.nodes.Delay(jnetwork, nodeName);
                elseif isa(node, 'Fork')
                    jnode = jline.lang.nodes.Fork(jnetwork, nodeName);
                elseif isa(node, 'Join')
                    % Note: Join requires Fork reference - simplified here
                    jnode = jline.lang.nodes.Join(jnetwork, nodeName, []);
                else
                    warning('NetworkXMLIO:UnknownNodeType', 'Unknown node type: %s', class(node));
                end
            end
        end
        
        function convertClasses(mnetwork, jnetwork)
            % CONVERTCLASSES(MNETWORK, JNETWORK)
            %
            % Converts job classes from MNetwork to JNetwork
            
            classes = mnetwork.getClasses();
            for i = 1:length(classes)
                jobClass = classes{i};
                className = jobClass.getName();
                priority = jobClass.getPriority();
                
                if isa(jobClass, 'OpenClass')
                    jclass = jline.lang.OpenClass(jnetwork, className, priority);
                elseif isa(jobClass, 'ClosedClass')
                    population = jobClass.getPopulation();
                    % Note: Reference station conversion needed
                    refStation = []; % Simplified
                    jclass = jline.lang.ClosedClass(jnetwork, className, population, refStation, priority);
                else
                    warning('NetworkXMLIO:UnknownClassType', 'Unknown class type: %s', class(jobClass));
                end
            end
        end
        
        function convertRoutingMatrices(mnetwork, jnetwork)
            % CONVERTROUTINGMATRICES(MNETWORK, JNETWORK)
            %
            % Converts routing matrices including class switching
            
            % Note: This is a placeholder for routing matrix conversion
            % Full implementation would need to:
            % 1. Extract routing probabilities P{i,j}(m,n) for all class combinations
            % 2. Handle class switching matrices from ClassSwitch nodes
            % 3. Convert to Java RoutingMatrix format
            
            warning('NetworkXMLIO:RoutingConversionTodo', ...
                   'Routing matrix conversion with class switching not fully implemented');
        end
        
        function convertServiceDistributions(mnetwork, jnetwork)
            % CONVERTSERVICEDISTRIBUTIONS(MNETWORK, JNETWORK)
            %
            % Converts service distributions from MNetwork to JNetwork
            
            % Note: This is a placeholder for service distribution conversion
            % Full implementation would need to:
            % 1. Extract service distributions for each node and class
            % 2. Convert MATLAB distribution objects to Java equivalents
            % 3. Handle arrival processes for Source nodes
            
            warning('NetworkXMLIO:ServiceConversionTodo', ...
                   'Service distribution conversion not fully implemented');
        end
        
        function result = yesNo(value)
            % YESNO(VALUE)
            %
            % Converts boolean to Yes/No string
            
            if value
                result = 'Yes';
            else
                result = 'No';
            end
        end
        
    end
end