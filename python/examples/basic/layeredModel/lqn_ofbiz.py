# %%
from line_solver import *
import numpy as np
import time
import os
GlobalConstants.set_verbose(VerboseLevel.SILENT)

# This example is temporarily disabled (matching MATLAB comment)
# Complex XML parsing may have incomplete support
# Clear variables (equivalent to MATLAB's clear solver AvgTable)
solver = None
AvgTable = None
# %%
print('This example illustrates the solution of a moderately large LQN.')

# MATLAB: cwd = fileparts(which(mfilename));
# MATLAB: model = LayeredNetwork.parseXML([cwd,filesep,'ofbizExample.xml']);
cwd = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd()
xml_path = os.path.join(cwd, 'lqn_ofbiz.xml')
#xml_path = os.path.join(cwd, 'lqn_ofbiz_fcfs.xml')

try:
    # Parse XML file like MATLAB: model = LayeredNetwork.parse_xml(filename)
    model = LayeredNetwork.parse_xml(xml_path)
    print(f'Successfully loaded model from {xml_path}')
except Exception as e:
    print(f'XML parsing failed: {e}')
    model = None
# %%
if model is None:
    print('Model could not be loaded, skipping solver execution')
    exit(0)

# MATLAB: options = LQNS.defaultOptions;
# MATLAB: options.keep = true; % uncomment to keep the intermediate XML files
#options = LQNS.default_options()
#options.keep = True

# Solve with LQNS (disabled)
# MATLAB: solver{1} = LQNS(model);
# MATLAB: AvgTable{1} = solver{1}.getAvgTable;
solver = {}
AvgTable = {}

#if LQNS.isAvailable():
#    try:
#        solver[1] = LQNS(model)
#        AvgTable[1] = solver[1].avg_table()
#        print('\nLQNS Results:')
#        print(AvgTable[1])
#    except Exception as e:
#        print(f'LQNS solver failed: {e}')
#else:
#    print('LQNS solver not available - skipping solver[1]')

# Solve with LN without initialization
# MATLAB: solver{2} = LN(model, @(x) NC(x,'verbose',false), lnoptions);
# MATLAB: AvgTable{2} = solver{2}.getAvgTable;
try:
    nc_options = NC.default_options()
    nc_options.verbose = False
    # Pin the LN seed, as the MATLAB twin does: SolverLN randomizes its
    # stochastic-layer seed base when options.seed is unset, which makes the
    # fixed point of this large model irreproducible run to run.
    lnoptions = LN.default_options()
    lnoptions.seed = 23000
    solver[2] = LN(model, lambda x: NC(x, nc_options), lnoptions)
    Tnoinit_start = time.time()
    AvgTable[2] = solver[2].avg_table()
    Tnoinit = time.time() - Tnoinit_start
    print('\nLN(NC) Results:')
    print(AvgTable[2])
    print(f'Tnoinit = {Tnoinit:.6f}')
except Exception as e:
    print(f'LN(NC) solver failed: {e}')
