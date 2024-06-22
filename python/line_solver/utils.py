import line_solver
from line_solver import GlobalConstants, VerboseLevel


import pandas as pd

def tget(table, *args):
    column = None
    results = None

    if len(args) == 1 or len(args) == 2:
        if args[0] in table.columns:
            column = args[0]
            remaining_args = args[1:]
        else:
            remaining_args = args[:]
    else:
        raise ValueError("Invalid number of arguments.")

    if len(remaining_args) == 1:
        if isinstance(remaining_args[0], str):  # default is Station
            value = remaining_args[0]
            results = table.loc[(table["Station"] == value)]
            if results.empty:
                results = table.loc[(table["JobClass"] == value)]
        elif isinstance(remaining_args[0], line_solver.lang.JobClass):
            jobclass = remaining_args[0].obj.getName()
            results = table.loc[(table["JobClass"] == jobclass)]
        else:
            station = remaining_args[0].obj.getName()
            results = table.loc[(table["Station"] == station)]
    elif len(remaining_args) == 2:
        if isinstance(remaining_args[0], str):  # default is Station
            station = remaining_args[0]
            jobclass = remaining_args[1]
            results = table.loc[(table["Station"] == station) & (table["JobClass"] == jobclass)]
        elif isinstance(remaining_args[1], line_solver.lang.JobClass):
            station = remaining_args[0].obj.getName()
            jobclass = remaining_args[1].obj.getName()
            results = table.loc[(table["Station"] == station) & (table["JobClass"] == jobclass)]

#    if not (GlobalConstants.getVerbose() == VerboseLevel.SILENT):
#        print(results)
    
    if column:
        if results is not None and not results.empty:
            if column in results.columns:
                return results[column].values[0]
            else:
                raise ValueError(f"Performance metric {column} is not supported.")
        else:
            return None
    else:
        return results
