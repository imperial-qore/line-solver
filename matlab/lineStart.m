%% LINE Initialization Script
% This script initializes the LINE environment by:
% - Setting global constants
% - Configuring default solver and tolerance parameters
% - Loading Java dependencies
% - Printing environment information

%% --- Global Variables ---
global LINEStdOut LINEVerbose LINEVersion LINEDoChecks LINEDummyMode
global LINECoarseTol LINEFineTol LINEImmediate LINEZero LINEMaxInt
global BuToolsVerbose BuToolsCheckInput BuToolsCheckPrecision
global LINELibraryAttributionShown LINEToolAckShown
global LINEDefaultLang

%% --- Setup Workspace Environment ---
cwd = fileparts(mfilename('fullpath'));
addpath(genpath(cwd)); % Add all subfolders to path
addpath(genpath([cwd,filesep,'..',filesep,'common',filesep])); % Add all subfolders to path

format compact

%% --- Clean up temporary files from past executions ---
try
    warning off
    lineClearWorkspace;
catch
    % Ignore errors if workspace folder doesn't exist yet
end
warning on backtrace

%% --- Default Configuration Values ---
LINE_VERSION     = '3.0.7';
LINE_STDOUT      = 1;        % 1 = console
LINE_VERBOSE     = VerboseLevel.STD;
LINE_DO_CHECKS   = true;
LINE_DUMMY_MODE  = false;
COARSE_TOL       = 1e-3;
FINE_TOL         = 1e-8;
ZERO_THRESHOLD   = 1e-14;
MAX_INT          = 2147483647; % same as Java Integer.MAX_VALUE

BUT_VERBOSE      = false;
BUT_CHECK_INPUT  = false;
BUT_PRECISION    = 1e-12;

%% --- Java Dependencies ---
% Get path to common folder
dev_dir = fileparts(mfilename('fullpath'));
root_dir = fileparts(dev_dir);
common_dir = fullfile(root_dir, 'common');

% Resolve the jline JAR path. The JLINE_JAR environment variable provides an
% explicit override for testing alternative builds;
% otherwise default to jline.jar in the common directory.
jlineJarPath = getenv('JLINE_JAR');
if isempty(jlineJarPath)
    jlineJarPath = fullfile(common_dir, 'jline.jar');
end
if ~exist(jlineJarPath, 'file')
    error('jline JAR not found at %s', jlineJarPath);
end

% ldes.jar is not loaded: its shaded jline copy would shadow jline.jar. MATLAB
% uses the LDES bundled in jline.jar; only python consumes ldes.jar.
if ~any(strcmp(javaclasspath('-dynamic'), jlineJarPath))
    javaaddpath(jlineJarPath);
end

import org.ejml.*; %#ok<SIMPT>

%% --- Assign Global Settings ---
LINEVersion     = LINE_VERSION;
LINEStdOut      = LINE_STDOUT;
LINEVerbose     = LINE_VERBOSE;
LINEDoChecks    = LINE_DO_CHECKS;
LINEDummyMode   = LINE_DUMMY_MODE;
LINECoarseTol   = COARSE_TOL;
LINEFineTol     = FINE_TOL;
LINEImmediate   = 1 / FINE_TOL;
LINEZero        = ZERO_THRESHOLD;
LINEMaxInt      = MAX_INT;

BuToolsVerbose        = BUT_VERBOSE;
BuToolsCheckInput     = BUT_CHECK_INPUT;
BuToolsCheckPrecision = BUT_PRECISION;

LINELibraryAttributionShown = false;
LINEToolAckShown = struct();

%% --- Startup Information ---
printLineStartup();

%% --- Local Utility Functions ---
function printLineStartup()
    global LINEStdOut LINEVersion LINEVerbose LINEDoChecks
    global LINECoarseTol LINEFineTol LINEZero LINEMaxInt

    fprintf(1, 'Starting LINE version %s: ', LINEVersion);

    if LINEStdOut == 1
        fprintf(1, 'StdOut=console, ');
    end

    switch LINEVerbose
        case VerboseLevel.STD
            fprintf('VerboseLevel=STD, ');
        case VerboseLevel.DEBUG
            fprintf('VerboseLevel=DEBUG, ');
        case VerboseLevel.DISABLED
            fprintf('VerboseLevel=DISABLED, ');
    end

    fprintf('DoChecks=%s, ', mat2str(LINEDoChecks));
    fprintf('CoarseTol=%.1e, FineTol=%.1e, Zero=%.1e, MaxInt=%d\n',...
        LINECoarseTol, LINEFineTol, LINEZero, LINEMaxInt);
    % Attribution is pull-based, as in Sage: the banner carries a pointer of
    % the same shape as Sage's 'Type "help()" for help.', and nothing is
    % printed during a solve.
    fprintf(1, ['Type solver.libraries() for third-party dependencies, ' ...
        'solver.citations() for references.\n']);
end
