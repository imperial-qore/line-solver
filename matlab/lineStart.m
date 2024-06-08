global LINEStdOut
global LINEVerbose
global LINEVersion
global LINEDoChecks
global LINEDummyMode % If true, getAvgTable returns empty results without running the solvers
global LINECoarseTol % Coarse tolerance for comparing averages, weights, states
global LINEFineTol % Fine tolerance eg for distribution comparisons
global LINEImmediate % Representation of infinite rate
global LINEZero % Generic representation of zero
global BuToolsVerbose
global BuToolsCheckInput
global BuToolsCheckPrecision

% keep this block here
cwd = fileparts(mfilename('fullpath'));
addpath(genpath(cwd));
format compact
warning ON BACKTRACE

% import java classes
javaaddpath(which('jline.jar'));
import org.ejml.*; %#ok<SIMPT>
javaaddpath(which('pfqn_nclib.jar'));
import DataStructures.*; %#ok<SIMPT>
import QueueingNet.*; %#ok<SIMPT>
import Utilities.*; %#ok<SIMPT>

% assign global constants
LINEStdOut = 1; % console
LINEVersion = '2.0.35';
LINEVerbose = VerboseLevel.STD;
LINEDoChecks = true;
LINEDummyMode = false; 
LINECoarseTol = 1e-3; 
LINEFineTol = 1e-8; 
LINEImmediate = 1/LINEFineTol;
LINEZero = 1e-14;

fprintf(1,'Starting LINE version %s: ',LINEVersion);
switch LINEStdOut
    case 1
    fprintf(1,'StdOut=console, ');
end 

switch LINEVerbose
    case VerboseLevel.STD
    fprintf(1,'VerboseLevel=STD, ');
    case VerboseLevel.DEBUG
    fprintf(1,'VerboseLevel=DEBUG, ');
    case VerboseLevel.DISABLED
    fprintf(1,'VerboseLevel=DISABLED, ');    
end

if LINEDoChecks
    fprintf(1,'DoChecks=true, ');
else
    fprintf(1,'DoChecks=false, ');
end

fprintf(1,'CoarseTol=%e, ',LINECoarseTol);
fprintf(1,'FineTol=%e, ',LINEFineTol);
fprintf(1,'Zero=%e\n',LINEZero);

BuToolsVerbose = false;
BuToolsCheckInput = false;
BuToolsCheckPrecision = 1e-12;
