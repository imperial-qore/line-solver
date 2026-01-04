"""
Utility functions for LINE queueing network analysis.

This module provides helper functions for working with LINE models
and results, including table manipulation, mathematical utilities,
and data processing functions.
"""

import pandas as pd
import numpy as np
from scipy.linalg import circulant
from .io.indexed_table import IndexedTable

__all__ = ['tget', 'circul', 'IndexedTable']


def tget(df, *args):
    """
    Extract specific rows/columns from LINE result tables.

    This function filters and selects data from pandas DataFrames containing
    LINE solver results based on station names, job class names, or other
    identifiers.

    Args:
        df (pandas.DataFrame): Input DataFrame with LINE results.
        *args: Variable arguments specifying filters (station names,
              job classes, or other identifiers).

    Returns:
        pandas.DataFrame: Filtered DataFrame with selected rows and columns.

    Examples:
        >>> results = solver.avg_table()
        >>> queue_results = tget(results, 'Queue')
        >>> class1_results = tget(results, 'Class1')
    """
    if not args:
        return df

    mask = pd.Series([True] * len(df), index=df.index)

    columns = df.columns.tolist()

    default_columns = ['Station', 'JobClass']

    for arg in args:
        if hasattr(arg, 'getName'):
            arg_value = str(arg.get_name())
        else:
            arg_value = str(arg)

        if arg_value in df.columns:
            columns = default_columns + [arg_value]
        else:
            mask = mask & df.apply(lambda row: row.astype(str).str.contains(arg_value).any(), axis=1)

    return df.loc[mask, columns].drop_duplicates()


def circul(c):
    """
    Generate a circulant matrix.

    Creates a circulant matrix where each row is a cyclic permutation of
    the previous row. For a scalar input, creates a circulant matrix of
    size c x c with specific pattern.

    Args:
        c: Either an integer (size of matrix) or array-like (first row).

    Returns:
        numpy.ndarray: The circulant matrix.

    Examples:
        >>> circul(3)  # Creates a 3x3 circulant matrix
        >>> circul([1, 2, 3])  # Creates circulant matrix with [1,2,3] as first row
    """
    if isinstance(c, (int, float)):
        n = int(c)
        # Create first row as [0, 1, 0, 0, ..., 0]
        first_row = np.zeros(n)
        if n > 1:
            first_row[1] = 1
        return circulant(first_row)
    else:
        c_arr = np.asarray(c).flatten()
        return circulant(c_arr)
