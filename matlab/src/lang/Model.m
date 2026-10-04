classdef Model < Copyable
    % Abstract parent class for all models in the LINE framework
    %
    % This class provides the basic structure and common functionality
    % for all LINE model types. It maintains model metadata such as
    % name, version, and attributes.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Hidden)
        attribute;      % Model attributes and metadata
        lineVersion;    % Version of LINE used to create the model
    end

    properties
        name;           % Name of the model
    end

    methods
        function self = Model(name)
            % Constructor: Creates a new Model instance
            %
            % Input:
            %   name - String name for the model
            %
            % Output:
            %   self - Model instance
            %
            % This constructor ensures LINE is properly initialized,
            % sets the model name, and records the LINE version.

            % Initialize LINE if not already done
            if isempty(GlobalConstants.Verbose)
                lineStart
            end
            
            % Set version and name
            lineVersion = strtrim(GlobalConstants.Version);
            self.setVersion(lineVersion);
            self.setName(name);
        end

        function out = getName(self)
            % Get the model name
            %
            % Output:
            %   out - String name of the model

            out = self.name;
        end

        function self = setName(self, name)
            % Set the model name
            %
            % Input:
            %   name - String name for the model
            %
            % Output:
            %   self - Updated Model instance
            
            self.name = name;
        end

        function v = getVersion(self)
            % Get the LINE version used to create this model
            %
            % Output:
            %   v - String version of LINE
            
            v = self.lineVersion;
        end

        function self = setVersion(self, version)
            % Set the LINE version for this model
            %
            % Input:
            %   version - String version of LINE
            %
            % Output:
            %   self - Updated Model instance

            self.lineVersion = version;
        end
    end

end
