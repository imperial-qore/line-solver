import line_solver
from line_solver import GlobalConstants, VerboseLevel


def tget(table, *argv):
    if len(argv) == 1:
        if isinstance(argv[0], str):  # default is station
            value = argv[0]
            results = table.loc[(table["Station"] == value)]
            if results.empty:
                results = table.loc[(table["JobClass"] == value)]
        elif isinstance(argv[0], line_solver.lang.JobClass):
            jobclass = argv[0].obj.getName()
            results = table.loc[(table["JobClass"] == jobclass)]
        else:
            station = argv[0].obj.getName()
            results = table.loc[(table["Station"] == station)]
    elif len(argv) == 2:
        if isinstance(argv[0], str):  # default is station
            station = argv[0]
            jobclass = argv[1]
            results = table.loc[(table["Station"] == station) & (table["JobClass"] == jobclass)]
        elif isinstance(argv[1], line_solver.lang.JobClass):
            station = argv[0].obj.getName()
            jobclass = argv[1].obj.getName()
            results = table.loc[(table["Station"] == station) & (table["JobClass"] == jobclass)]
    if not (GlobalConstants.getVerbose() == VerboseLevel.SILENT):
        print(results)
    return results
