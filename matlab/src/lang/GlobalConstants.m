classdef GlobalConstants
    % GlobalConstants System-wide constants and configuration parameters
    %
    % GlobalConstants provides centralized access to global constants, tolerances,
    % and configuration parameters used throughout the LINE framework. It manages
    % numerical tolerances, verbosity levels, version information, and other
    % system-wide settings through static methods and global variables.
    %
    % @brief Centralized global constants and configuration management
    %
    % Key characteristics:
    % - Centralized constant management
    % - Numerical tolerance configuration
    % - System-wide parameter access
    % - Version and build information
    % - Debugging and verbosity controls
    %
    % Global constants include:
    % - Numerical tolerances (FineTol, CoarseTol)
    % - Special values (Zero, MaxInt, Immediate)
    % - System configuration (Verbose, DummyMode)
    % - Version information and build details
    % - Output stream redirection (StdOut)
    %
    % GlobalConstants is used for:
    % - Consistent numerical precision across LINE
    % - Global configuration management
    % - Debugging and verbosity control
    % - Version compatibility checking
    % - System-wide parameter standardization
    %
    % Example:
    % @code
    % if abs(value) < GlobalConstants.FineTol()
    %     % Handle near-zero values
    % end
    % if GlobalConstants.Verbose() > 1
    %     fprintf('Debug information\n');
    % end
    % @endcode
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods (Static)
        function can=DummyMode()
            global LINEDummyMode
            can = LINEDummyMode;
        end

        function stdo=StdOut()
            global LINEStdOut
            stdo = LINEStdOut;
        end

        function tol=Immediate()
            global LINEImmediate
            tol = LINEImmediate;
        end

        function tol=MaxInt()
            global LINEMaxInt
            tol = LINEMaxInt;
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

        function setChecks(val)
            global LINEChecks
            LINEChecks = val;
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

        function setMaxInt(val)
            global LINEMaxInt
            LINEMaxInt = val;
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

        function verbose = getVerbose()
            global LINEVerbose
            verbose = LINEVerbose;
        end

        function setVersion(val)
            global LINEVersion
            LINEVersion = val;
        end

        function bool = isLibraryAttributionShown()
            global LINELibraryAttributionShown
            bool = LINELibraryAttributionShown;
        end

        function setLibraryAttributionShown(val)
            global LINELibraryAttributionShown
            LINELibraryAttributionShown = val;
        end
    end
end
