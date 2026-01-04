
"""
Polling system analysis functions.

This module provides analytical methods for polling systems where
a single server visits multiple queues according to a polling protocol.
Supports various polling disciplines including exhaustive, gated,
and k-limited service.

These functions analyze systems with switchover times and
different service disciplines at each queue.
"""

import jpype
import numpy as np
from line_solver import jlineMatrixToArray, jlineMatrixFromArray


def polling_qsys_1limited(arv_maps, svc_maps, switch_maps):
    """
    Analyze 1-limited polling system.
    
    Models a polling system where the server serves at most one customer
    per visit to each queue before moving to the next queue.
    
    Args:
        arv_maps: List of arrival MAPs for each queue.
        svc_maps: List of service MAPs for each queue.
        switch_maps: List of switchover time MAPs between queues.
        
    Returns:
        dict: Performance measures for the 1-limited polling system.
    """
    def convert_to_matrix_cell_array(maps_list):
        java_array = jpype.JArray(jpype.JPackage('jline').util.matrix.MatrixCell)(len(maps_list))
        for i, map_item in enumerate(maps_list):
            if isinstance(map_item, tuple):
                D0, D1 = map_item
                java_cell = jpype.JPackage('jline').util.matrix.MatrixCell(2)
                java_cell.add(jlineMatrixFromArray(D0))
                java_cell.add(jlineMatrixFromArray(D1))
                java_array[i] = java_cell
            else:
                java_array[i] = map_item
        return java_array

    java_arv_maps = convert_to_matrix_cell_array(arv_maps)
    java_svc_maps = convert_to_matrix_cell_array(svc_maps)
    java_switch_maps = convert_to_matrix_cell_array(switch_maps)

    result = jpype.JPackage('jline').api.polling.Polling_qsys_1limitedKt.polling_qsys_1limited(
        java_arv_maps, java_svc_maps, java_switch_maps
    )

    return np.array(result)


def polling_qsys_exhaustive(arv_maps, svc_maps, switch_maps):
    """
    Analyze exhaustive service polling system.

    Args:
        arv_maps: List of arrival MAPs for each queue.
        svc_maps: List of service MAPs for each queue.
        switch_maps: List of switchover time MAPs between queues.

    Returns:
        numpy.ndarray: Performance measures for exhaustive polling system.
    """
    def convert_to_matrix_cell_array(maps_list):
        java_array = jpype.JArray(jpype.JPackage('jline').util.matrix.MatrixCell)(len(maps_list))
        for i, map_item in enumerate(maps_list):
            if isinstance(map_item, tuple):
                D0, D1 = map_item
                java_cell = jpype.JPackage('jline').util.matrix.MatrixCell(2)
                java_cell.add(jlineMatrixFromArray(D0))
                java_cell.add(jlineMatrixFromArray(D1))
                java_array[i] = java_cell
            else:
                java_array[i] = map_item
        return java_array

    java_arv_maps = convert_to_matrix_cell_array(arv_maps)
    java_svc_maps = convert_to_matrix_cell_array(svc_maps)
    java_switch_maps = convert_to_matrix_cell_array(switch_maps)

    result = jpype.JPackage('jline').api.polling.Polling_qsys_exhaustiveKt.polling_qsys_exhaustive(
        java_arv_maps, java_svc_maps, java_switch_maps
    )

    return np.array(result)


def polling_qsys_gated(arv_maps, svc_maps, switch_maps):
    """
    Analyze gated service polling system.

    Args:
        arv_maps: List of arrival MAPs for each queue.
        svc_maps: List of service MAPs for each queue.
        switch_maps: List of switchover time MAPs between queues.

    Returns:
        numpy.ndarray: Performance measures for gated polling system.
    """
    def convert_to_matrix_cell_array(maps_list):
        java_array = jpype.JArray(jpype.JPackage('jline').util.matrix.MatrixCell)(len(maps_list))
        for i, map_item in enumerate(maps_list):
            if isinstance(map_item, tuple):
                D0, D1 = map_item
                java_cell = jpype.JPackage('jline').util.matrix.MatrixCell(2)
                java_cell.add(jlineMatrixFromArray(D0))
                java_cell.add(jlineMatrixFromArray(D1))
                java_array[i] = java_cell
            else:
                java_array[i] = map_item
        return java_array

    java_arv_maps = convert_to_matrix_cell_array(arv_maps)
    java_svc_maps = convert_to_matrix_cell_array(svc_maps)
    java_switch_maps = convert_to_matrix_cell_array(switch_maps)

    result = jpype.JPackage('jline').api.polling.Polling_qsys_gatedKt.polling_qsys_gated(
        java_arv_maps, java_svc_maps, java_switch_maps
    )

    return np.array(result)