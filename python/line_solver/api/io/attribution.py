"""
Once-per-session attribution of the third-party libraries a solver uses.

This is the counterpart of MATLAB's @NetworkSolver/showLibraryAttribution.m and
of jline.solvers.NetworkSolver.showLibraryAttribution: LINE bundles ports of
BUTools, Q-MAM, MAMSolver, rmf_tool and others under lib/thirdparty, and names
them once per session so a user knows whose algorithms produced the numbers.
It is distinct from line_ack, which credits an external TOOL a wrapper solver
shells out to (JMT, LQNS, qnsolver).
"""

from ...constants import GlobalConstants


def libraries_used(solver_name, sn=None, options=None):
    """
    Third-party libraries a solver will use, as MATLAB's per-solver static
    getLibrariesUsed reports them.

    Args:
        solver_name: short solver name, e.g. 'MAM' or 'FLD'
        sn: NetworkStruct of the model, or None
        options: solver options, or None

    Returns:
        list of library names (possibly empty)
    """
    libs = []
    method = ''
    if options is not None:
        method = str(getattr(options, 'method', '') or '')
    name = (solver_name or '').upper()

    if name == 'MAM':
        if method in ('default', 'dec.source', 'dec.mmap', 'dec.poisson', 'dec.source.mmap'):
            libs.append('MAMSolver')
        if method in ('mna', 'inap', 'inapplus', 'inapinf'):
            libs.append('Q-MAM')
        if sn is not None and getattr(sn, 'proc', None):
            libs.append('BUTools')
    elif name in ('FLD', 'FLUID'):
        if method in ('rmf', 'fluid.rmf'):
            libs.append('rmf_tool')
    return libs


def show_library_attribution(solver_name, sn=None, options=None):
    """
    Print the attribution once per session, honouring the verbosity setting.

    Args:
        solver_name: short solver name, e.g. 'MAM'
        sn: NetworkStruct of the model, or None
        options: solver options, or None
    """
    from ...constants import VerboseLevel
    verbose = getattr(options, 'verbose', None) if options is not None else None
    if verbose is not None and verbose == VerboseLevel.SILENT:
        return
    if GlobalConstants.isLibraryAttributionShown():
        return
    libs = libraries_used(solver_name, sn, options)
    if not libs:
        return
    print('The solver will leverage %s.' % ', '.join(libs))
    GlobalConstants.setLibraryAttributionShown(True)


def solver_libraries(solver):
    """
    Third-party libraries a solver instance will use, without printing.

    Attribution in LINE is pull-based, in the spirit of Sage's
    sage.misc.citation.get_systems: nothing is written to the console during a
    solve, and the user asks for the list when citing.

    Args:
        solver: a solver instance exposing .options and, optionally, a struct

    Returns:
        list of library names, possibly empty
    """
    name = type(solver).__name__.replace('Solver', '')
    sn = getattr(solver, 'sn', None)
    if sn is None:
        sn = getattr(solver, '_sn', None)
    return libraries_used(name, sn, getattr(solver, 'options', None))
