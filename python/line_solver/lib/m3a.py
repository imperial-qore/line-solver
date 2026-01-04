"""
3rd Moment Approximation (M3A) utilities.

Provides functions for fitting and using 3rd-moment approximations:
- Auto-fitting of M3A parameters
- Manual parameter specification
- Compression and fitting of arrival processes
- MMAP (Marked Markovian Arrival Process) fitting and compression

References:
[1] A. Sansottera, G. Casale, P. Cremonesi. Fitting Second-Order Acyclic
    Marked Markovian Arrival Processes. IEEE/IFIP DSN 2013.
[2] G. Casale, A. Sansottera, P. Cremonesi. Compact Markov-Modulated
    Models for Multiclass Trace Fitting. European Journal of Operations
    Research, 2016.
"""

import numpy as np
from line_solver import native_to_array


def _list_to_matrix_cell(mmap_list):
    """Convert a Python list of numpy arrays to MatrixCell."""
    cell = MatrixCell(len(mmap_list))
    for i, mat in enumerate(mmap_list):
        cell.set(i, native_to_array(np.asarray(mat)))
    return cell


def _matrix_cell_to_list(cell):
    """Convert MatrixCell to a Python list of numpy arrays."""
    result = []
    for i in range(cell.size()):
        result.append(native_to_array(cell.get(i)))
    return result


def m3afit_init(S, C):
    """
    Prepare multiclass trace for M3A fitting.

    Args:
        S: Inter-arrival times (array-like)
        C: Class number for each arrival (array-like of integers)

    Returns:
        tuple: (S, C, num_classes) ready for m3afit_auto
    """
    S = np.asarray(S).flatten()
    C = np.asarray(C).flatten().astype(int)
    num_classes = len(np.unique(C))
    return (S, C, num_classes)


def m3afit_auto(mtrace, num_states=2, method=1, timescale=None, timescale_asy=None):
    """
    Automatic fitting of trace into a Marked Markovian Arrival Process.

    Based on the M3A toolbox, this function selects the appropriate fitting
    algorithm based on the number of classes, requested states, and fitting method.

    Args:
        mtrace: Tuple (S, C, num_classes) returned by m3afit_init, or just (S, C)
        num_states: Number of states for the fitted MMAP (default: 2)
        method: Fitting method (0 = inter-arrival, 1 = counting process)
        timescale: Finite time scale for counting process (auto-computed if None)
        timescale_asy: Near-infinite time scale (auto-computed if None)

    Returns:
        list: Fitted MMAP as [D0, D1, D11, D12, ...] or None if fitting fails

    Example:
        >>> S = [1.0, 0.5, 1.2, 0.8, 0.3]  # inter-arrival times
        >>> C = [0, 1, 0, 1, 0]  # class labels
        >>> mtrace = m3afit_init(S, C)
        >>> mmap = m3afit_auto(mtrace, num_states=2, method=1)
    """
    try:

        # Extract S and C from mtrace
        if len(mtrace) == 3:
            S, C, _ = mtrace
        else:
            S, C = mtrace

        S = np.asarray(S).flatten()
        C = np.asarray(C).flatten().astype(int)

        # Convert to arrays

        # Call function
        result_cell = m3afit_auto_fn(S_arr, C_arr, int(num_states), int(method))

        if result_cell is None:
            return None

        return _matrix_cell_to_list(result_cell)
    except Exception as e:
        raise RuntimeError(f"M3A auto fitting failed: {str(e)}")


def m3afit_compress(mmap, method=0, num_states=2, verbose=False):
    """
    Compress a Marked Markovian Arrival Process using M3A approximation.

    Takes an arbitrary-order MMAP and produces a compressed second-order
    acyclic MMAP that approximates the original process.

    The compression preserves key statistical characteristics:
    - First three moments of inter-arrival times
    - Autocorrelation structure (via gamma decay rate)
    - Class probabilities
    - Forward and backward moments

    Args:
        mmap: Input MMAP as list [D0, D1, D11, D12, ...] where:
              - D0: Transition matrix without arrivals
              - D1: Aggregate arrival matrix (sum of D1c)
              - D1c: Class c arrival matrix for c = 0, 1, ...
        method: Compression method (0 = 2-state AMAP compression)
        num_states: Number of states in compressed representation
        verbose: Print progress information

    Returns:
        list: Compressed MMAP as [D0, D1, D11, D12, ...]

    Example:
        >>> D0 = np.array([[-1.0, 0.2], [0.3, -0.8]])
        >>> D1 = D11 + D12  # sum of class matrices
        >>> D11 = np.array([[0.3, 0.1], [0.1, 0.2]])
        >>> D12 = np.array([[0.2, 0.1], [0.1, 0.2]])
        >>> mmap = [D0, D1, D11, D12]
        >>> compressed = m3afit_compress(mmap)
    """
    try:

        # Convert Python list to MatrixCell
        mmap_cell = _list_to_matrix_cell(mmap)

        # Call function
        result_cell = m3afit_compress_fn(mmap_cell)

        # Convert result back to Python list
        return _matrix_cell_to_list(result_cell)
    except Exception as e:
        raise RuntimeError(f"M3A compression failed: {str(e)}")


def mamap2m_fit_gamma_fb(M1, M2, M3, GAMMA, P, F, B):
    """
    Fit a MAMAP(2,m) matching moments, autocorrelation, and class characteristics.

    Computes a second-order MAMAP[m] fitting the given ordinary moments,
    autocorrelation decay rate, class probabilities, forward moments,
    and backward moments.

    Args:
        M1: First moment of inter-arrival times
        M2: Second moment of inter-arrival times
        M3: Third moment of inter-arrival times
        GAMMA: Autocorrelation decay rate
        P: Class probabilities (array of length m)
        F: First-order forward moments (array of length m)
        B: First-order backward moments (array of length m)

    Returns:
        list: Fitted MAMAP as [D0, D1, D11, D12, ...]
    """
    try:


        result_cell = fit_fn(float(M1), float(M2), float(M3), float(GAMMA), P_arr, F_arr, B_arr)

        return _matrix_cell_to_list(result_cell)
    except Exception as e:
        raise RuntimeError(f"MAMAP2M gamma FB fitting failed: {str(e)}")


def mamap2m_fit_gamma_fb_mmap(mmap):
    """
    Fit a second-order MAMAP from an existing MMAP.

    Extracts characteristics from the input MMAP and fits a second-order
    acyclic MMAP that approximates the original process.

    Args:
        mmap: Input MMAP as list [D0, D1, D11, D12, ...]

    Returns:
        list: Fitted second-order MAMAP as [D0, D1, D11, D12, ...]
    """
    try:

        mmap_cell = _list_to_matrix_cell(mmap)
        result_cell = fit_fn(mmap_cell)

        return _matrix_cell_to_list(result_cell)
    except Exception as e:
        raise RuntimeError(f"MAMAP2M gamma FB MMAP fitting failed: {str(e)}")


def mamap22_fit_gamma_bs(M1, M2, M3, GAMMA, P, B, S):
    """
    Fit a MAMAP(2,2) using backward moments and class transitions.

    Performs approximate fitting of a MMAP given the underlying MAP,
    the class probabilities, the backward moments, and the one-step
    class transition probabilities.

    Args:
        M1, M2, M3: Moments of inter-arrival times
        GAMMA: Autocorrelation decay rate
        P: Class probabilities (2-element array)
        B: First-order backward moments (2-element array)
        S: One-step class transition probabilities (2x2 matrix)

    Returns:
        list: Fitted MAMAP(2,2) as [D0, D1, D11, D12]
    """
    try:

        P_mat = native_to_array(np.asarray(P).reshape(1, -1))
        B_mat = native_to_array(np.asarray(B).reshape(-1, 1))
        S_mat = native_to_array(np.asarray(S))

        result_cell = fit_fn(float(M1), float(M2), float(M3), float(GAMMA), P_mat, B_mat, S_mat)

        return _matrix_cell_to_list(result_cell)
    except Exception as e:
        raise RuntimeError(f"MAMAP22 gamma BS fitting failed: {str(e)}")


def mamap22_fit_gamma_fs(M1, M2, M3, GAMMA, P, F, S):
    """
    Fit a MAMAP(2,2) using forward moments and class transitions.

    Performs approximate fitting of a MMAP given the underlying MAP,
    the class probabilities, the forward moments, and the one-step
    class transition probabilities.

    Args:
        M1, M2, M3: Moments of inter-arrival times
        GAMMA: Autocorrelation decay rate
        P: Class probabilities (2-element array)
        F: First-order forward moments (2-element array)
        S: One-step class transition probabilities (2x2 matrix)

    Returns:
        list: Fitted MAMAP(2,2) as [D0, D1, D11, D12]
    """
    try:

        P_mat = native_to_array(np.asarray(P).reshape(1, -1))
        F_mat = native_to_array(np.asarray(F).reshape(-1, 1))
        S_mat = native_to_array(np.asarray(S))

        result_cell = fit_fn(float(M1), float(M2), float(M3), float(GAMMA), P_mat, F_mat, S_mat)

        return _matrix_cell_to_list(result_cell)
    except Exception as e:
        raise RuntimeError(f"MAMAP22 gamma FS fitting failed: {str(e)}")


def mamap22_fit_gamma_bs_mmap(mmap):
    """
    Fit a MAMAP(2,2) from an existing MMAP using backward moments.

    Performs approximate fitting of an MMAP[2], yielding a second-order
    acyclic MMAP[2] fitting the backward moments and the class transition
    probabilities.

    Args:
        mmap: Input MMAP as list [D0, D1, D11, D12]

    Returns:
        list: Fitted MAMAP(2,2) as [D0, D1, D11, D12]
    """
    try:

        mmap_cell = _list_to_matrix_cell(mmap)
        result_cell = fit_fn(mmap_cell)

        return _matrix_cell_to_list(result_cell)
    except Exception as e:
        raise RuntimeError(f"MAMAP22 gamma BS MMAP fitting failed: {str(e)}")


def mamap22_fit_gamma_fs_mmap(mmap):
    """
    Fit a MAMAP(2,2) from an existing MMAP using forward moments.

    Performs approximate fitting of an MMAP[2], yielding a second-order
    acyclic MMAP[2] fitting the forward moments and the class transition
    probabilities.

    Args:
        mmap: Input MMAP as list [D0, D1, D11, D12]

    Returns:
        list: Fitted MAMAP(2,2) as [D0, D1, D11, D12]
    """
    try:

        mmap_cell = _list_to_matrix_cell(mmap)
        result_cell = fit_fn(mmap_cell)

        return _matrix_cell_to_list(result_cell)
    except Exception as e:
        raise RuntimeError(f"MAMAP22 gamma FS MMAP fitting failed: {str(e)}")


def amap2_fit_gamma(M1, M2, M3, GAMMA):
    """
    Fit an AMAP(2) (Acyclic MAP of order 2) from moments and autocorrelation.

    Args:
        M1: First moment of inter-arrival times
        M2: Second moment of inter-arrival times
        M3: Third moment of inter-arrival times
        GAMMA: Autocorrelation decay rate

    Returns:
        tuple: (best_map, all_maps) where:
            - best_map: Best-fit AMAP(2) as [D0, D1]
            - all_maps: List of all valid AMAP(2) representations
    """
    try:

        result = fit_fn(float(M1), float(M2), float(M3), float(GAMMA))

        # Result is a Pair<MatrixCell, List<MatrixCell>>
        best_map = _matrix_cell_to_list(result.getFirst())

        all_maps_list = result.getSecond()
        all_maps = []
        for i in range(all_maps_list.size()):
            all_maps.append(_matrix_cell_to_list(all_maps_list.get(i)))

        return best_map, all_maps
    except Exception as e:
        raise RuntimeError(f"AMAP2 gamma fitting failed: {str(e)}")


def lib_m3a_fit_auto(moments):
    """
    Automatically fit M3A parameters from moments.

    Args:
        moments (array-like): First 3 moments (mean, second, third)

    Returns:
        dict: M3A parameters
    """
    moments = np.asarray(moments).flatten()

    try:
        result = M3aFit.fitAuto(native_to_array(moments))

        return {
            'parameters': native_to_array(result),
            'method': 'auto'
        }
    except Exception as e:
        raise RuntimeError(f"M3A auto fitting failed: {str(e)}")


def lib_m3a_fit_init(moments):
    """
    Initialize M3A fitting process.

    Args:
        moments (array-like): First 3 moments

    Returns:
        dict: Initial M3A parameters
    """
    moments = np.asarray(moments).flatten()

    try:
        result = M3aFit.fitInit(native_to_array(moments))

        return {
            'parameters': native_to_array(result),
            'method': 'init'
        }
    except Exception as e:
        raise RuntimeError(f"M3A init fitting failed: {str(e)}")


def lib_m3a_compress(arrival_process):
    """
    Compress arrival process using M3A approximation.

    Args:
        arrival_process (array-like): Original arrival process representation

    Returns:
        dict: Compressed M3A representation with parameters
    """
    arrival_process = np.asarray(arrival_process)
    arrival_arr = native_to_array(arrival_process)

    try:

        result = M3aCompressor.compressHyperExponential(arrival_arr)

        return {
            'compressed': native_to_array(result),
            'method': 'hyperexponential'
        }
    except Exception as e:
        raise RuntimeError(f"M3A compression failed: {str(e)}")


def lib_m3a_summary(moments):
    """
    Provide summary of M3A fitting options for given moments.

    Args:
        moments (array-like): First 3 moments

    Returns:
        dict: Summary of M3A characteristics
    """
    moments = np.asarray(moments).flatten()

    m1 = moments[0] if len(moments) > 0 else 0
    m2 = moments[1] if len(moments) > 1 else 0
    m3 = moments[2] if len(moments) > 2 else 0

    scv = (m2 / (m1**2)) - 1 if m1 > 0 else 0
    skew = (m3 - 3*m2*m1 + 2*m1**3) / ((m2 - m1**2)**(3/2)) if (m2 - m1**2) > 0 else 0

    return {
        'mean': m1,
        'scv': scv,
        'skew': skew,
        'fitting_notes': f'SCV={scv:.3f}, Skew={skew:.3f}'
    }
