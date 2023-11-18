classdef GlobalConstants
    % Distribution for global constants
    %
    % Copyright (c) 2012-2023, Imperial College London
    % All rights reserved.

    methods (Static)
        function can=DummyMode()
            global LINEDummyMode
            can = LINEDummyMode;
        end

        function do=DoChecks()
            global LINEDoChecks
            do = LINEDoChecks;
        end

        function stdo=StdOut()
            global LINEStdOut
            stdo = LINEStdOut;
        end

        function tol=Immediate()
            global LINEImmediate
            tol = LINEImmediate;
        end

        function tol=Zero()
            global LINEZero
            tol = LINEZero;
        end        

        function tol=CoarseTol()
            global LINECoarseTol
            tol = LINECoarseTol;
        end

        function tol=FineTol()
            global LINEFineTol
            tol = LINEFineTol;
        end

        function verbose=Verbose()
            global LINEVerbose
            verbose = LINEVerbose;
        end

        function ver=Version()
            global LINEVersion
            ver = LINEVersion;
        end

        function setDoChecks(val)
            global LINEDoChecks
            LINEDoChecks = val;
        end

        function setStdOut(val)
            global LINEStdOut
            LINEStdOut = val;
        end

        function setImmediate(val)
            global LINEImmediate
            LINEImmediate = val;
        end

        function setDummyMode(val)
            global LINEDummyMode
            LINEDummyMode = val;
        end

        function setZero(val)
            global LINEZero
            LINEZero = val;
        end        

        function setCoarseTol(val)
            global LINECoarseTol
            LINECoarseTol = val;
        end

        function setFineTol(val)
            global LINEFineTol
            LINEFineTol = val;
        end

        function setVerbose(val)
            global LINEVerbose
            LINEVerbose = val;
        end

        function setVersion(val)
            global LINEVersion
            LINEVersion = val;
        end
    end
end
