# %%
from line_solver import *
import os

GlobalConstants.set_verbose(VerboseLevel.SILENT)

print('This example illustrates the initialization of LN using the output of LQNS.')

# Load model from XML file (in same directory as notebook)
model = LayeredNetwork.parse_xml('lqn_serial.xml')
# %%
# Set options for LQNS solver
options = LQNS.default_options()
options.keep = True  # Keep intermediate XML files

# Solve the model using LQNS if available
if LQNS.isAvailable():
    lqns_solver = LQNS(model, options)
    try:
        avg_table_lqns = lqns_solver.avg_table()
        print('LQNS Solver Results:')
        print(avg_table_lqns)
    except Exception as e:
        print(f'LQNS solver error: {e}')
else:
    print('LQNS solver not available - skipping')
# %%
# Solve with LN without initialization
print('\nSolve with LN without initialization:')
ln_solver = LN(model, lambda x: MVA(x))
import time
start_time = time.time()
try:
    avg_table = ln_solver.avg_table()
    end_time = time.time()
    time_elapsed = end_time - start_time
    
    print('AvgTable:')
    print(avg_table)
    print(f'Time elapsed: {time_elapsed:.3f}s')
except Exception as e:
    print(f'LN solver error: {e}')
# %%
# Obtain the CDF of response times
print('\nWe now obtain the CDF of response times:')
try:
    ensemble = model.ensemble
    if ensemble is not None and len(ensemble) >= 3:
        fluid_solver = FLD(ensemble[2])  # ensemble[2] in Python (0-indexed)
        rd = fluid_solver.cdf_resp_t()
        print('RD (CDF of response times):')
        print(rd)
    else:
        print('Model ensemble not available or insufficient layers')
except Exception as e:
    print(f'CDF computation not available: {e}')
    print('Fluid solver and ensemble access not yet fully implemented')