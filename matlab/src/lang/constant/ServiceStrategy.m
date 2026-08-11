classdef (Sealed) ServiceStrategy
    % Enumeration of service strategies
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties (Constant)
        LI = 1; % LoadIndependent
        LD = 2; % LoadDependent
        CD = 3; % ClassDependent
        SD = 4; % StateDependent
        JD = 5; % JointDependent
    end
    
    methods (Static)
        function txt = toText(strategy)
            switch strategy
                case ServiceStrategy.LI
                    txt = 'LoadIndependent';
                case ServiceStrategy.LD
                    txt = 'LoadDependent';
                case ServiceStrategy.CD
                    txt = 'ClassDependent';
                case ServiceStrategy.SD
                    txt = 'StateDependent';
                case ServiceStrategy.JD
                    txt = 'JointDependent';
                otherwise
                    line_error(mfilename, 'Unrecognized service strategy.');
            end
        end
    end

    methods (Access = private)
        function out = ServiceStrategy
            % Prevent instantiation
        end
    end
end
