classdef (Sealed) NodeType
    % Enumeration of node types.
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    properties (Constant)
        Region = 10;
        Transition = 9;
        Place = 8;
        Fork = 7;
        Router = 6;
        Cache = 5;
        Logger = 4;
        ClassSwitch = 3;
        Delay = 2;
        Source = 1;
        Queue = 0;
        Sink = -1;
        Join = -2;
    end

    methods (Static)

        function txt = toText(nodetype)
            switch nodetype
                case NodeType.Region
                    txt = 'Region';
                case NodeType.Transition
                    txt = 'Transition';
                case NodeType.Place
                    txt = 'Place';
                case NodeType.Fork
                    txt = 'Fork';
                case NodeType.Router
                    txt = 'Router';
                case NodeType.Cache
                    txt = 'Cache';
                case NodeType.Logger
                    txt = 'Logger';
                case NodeType.ClassSwitch
                    txt = 'ClassSwitch';
                case NodeType.Delay
                    txt = 'Delay';
                case NodeType.Source
                    txt = 'Source';
                case NodeType.Queue
                    txt = 'Queue';
                case NodeType.Join
                    txt = 'Join';
                case NodeType.Sink
                    txt = 'Sink';
                otherwise
                    line_error(mfilename, 'Unrecognized node type.');
            end
        end
    end

end
