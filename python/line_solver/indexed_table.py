"""
IndexedTable: Generic enhanced table wrapper with object-based filtering.

This module provides a pandas DataFrame wrapper that enables filtering using
Station, Node, JobClass, and Chain objects, while maintaining full compatibility
with standard pandas operations.

Supports all NetworkSolver result tables:
- AvgTable (Station + JobClass)
- AvgNodeTable (Node + JobClass)
- AvgChainTable (Station + Chain)
- AvgNodeChainTable (Node + Chain)
- AvgSysTable (Chain only)

Usage:
    from line_solver import *

    # Create and solve model
    model = Network('Example')
    # ... build model ...
    solver = SolverJMT(model)
    solver.runAnalyzer()

    # Get result table and wrap it
    avg_table_native = solver.avg_table()
    avg_table = IndexedTable(avg_table_native)

    # Four equivalent syntaxes for filtering
    result = avg_table[queue, jobclass]              # Direct indexing
    result = avg_table.filterBy(queue, jobclass)     # filterBy method
    result = avg_table.get(queue, jobclass)          # get alias
    result = avg_table.tget(queue, jobclass)         # tget alias

    # Access underlying DataFrame
    df = avg_table.data
"""

import pandas as pd
import numpy as np


def _matlab_format_float(x, precision=5):
    """
    Format a float in MATLAB style - uses fixed-point for moderate values,
    scientific notation for very small/large values.

    MATLAB uses 'g' style formatting that:
    - Uses fixed-point for values that can be represented compactly
    - Uses scientific notation for very small (<1e-4) or very large (>=1e5) values
    - Removes trailing zeros
    """
    if pd.isna(x):
        return 'NaN'
    if not isinstance(x, (int, float, np.integer, np.floating)):
        return str(x)

    # Handle special cases
    if x == 0:
        return '0'
    if np.isinf(x):
        return 'Inf' if x > 0 else '-Inf'

    abs_x = abs(x)

    # MATLAB uses scientific notation for very small or very large values
    # Use fixed-point if 1e-4 <= |x| < 1e5, otherwise scientific
    if abs_x >= 1e-4 and abs_x < 1e5:
        # Use fixed-point, try to show up to 'precision' significant figures
        # Format with enough decimal places, then strip trailing zeros
        formatted = f'{x:.{precision}g}'
    else:
        # Use scientific notation
        formatted = f'{x:.{precision-1}e}'

    return formatted


class IndexedTable:
    """
    Enhanced pandas DataFrame wrapper with object-based filtering.

    Wraps a MATLAB table to enable filtering using Station, Node, JobClass,
    and/or Chain objects, while maintaining full backward compatibility with
    standard pandas operations.

    Attributes:
        data (pd.DataFrame): The underlying pandas DataFrame
    """

    def __init__(self, dataframe):
        """
        Initialize IndexedTable wrapper.

        Args:
            dataframe: A pandas DataFrame (typically from solver.avgTable(), etc.)

        Raises:
            TypeError: If input is not a pandas DataFrame
        """
        if not isinstance(dataframe, pd.DataFrame):
            raise TypeError('IndexedTable requires a pandas DataFrame as input')
        self.data = dataframe.copy()

    def __getattr__(self, name):
        """
        Delegate attribute access to the underlying DataFrame.

        Columns are accessible as attributes (table.QLen, table.Util, ...) and
        any other DataFrame attribute/method (iterrows, to_string, values, ...)
        is forwarded so IndexedTable is a drop-in for standard pandas usage.
        """
        if name == 'data':
            raise AttributeError(name)
        if name in self.data.columns:
            return self.data[name]
        return getattr(self.data, name)

    def __getitem__(self, key):
        """
        Support direct indexing with objects: table[queue, jobclass]

        Args:
            key: Either a single object/tuple of objects for filtering,
                 or standard pandas indexing (int, slice, list, etc.)

        Returns:
            pd.DataFrame: Filtered result or standard pandas indexing result
        """
        # Handle tuple of objects for dual filtering
        if isinstance(key, tuple):
            if len(key) == 2 and self._is_object_arg(key[0]):
                # Object-based dual filtering
                return self.filterBy(key[0], key[1])
            else:
                # Standard pandas tuple indexing
                return self.data[key]

        # Handle single object filtering
        elif self._is_object_arg(key):
            return self.filterBy(key)

        # Standard pandas indexing
        else:
            return self.data[key]

    def filterBy(self, *args):
        """
        Filter table by Station/Node and/or JobClass/Chain objects.

        Intelligently filters based on table structure and object types.

        Args:
            *args: 1 or 2 arguments - Station/Node/JobClass/Chain objects

        Returns:
            pd.DataFrame: Filtered DataFrame

        Raises:
            ValueError: If invalid argument types or counts provided
        """
        if len(args) == 1:
            # Single object argument
            obj = args[0]
            if self._is_station_or_node(obj):
                return self._filter_by_first_dimension(obj)
            elif self._is_class_or_chain(obj):
                return self._filter_by_second_dimension(obj)
            else:
                # If object type can't be determined but has getName/name,
                # return empty DataFrame (graceful handling for missing columns)
                if hasattr(obj, 'getName') or hasattr(obj, 'name'):
                    return self.data.iloc[0:0].copy()  # Empty DataFrame with same structure
                raise ValueError(
                    'Invalid argument type. Expected Station, Node, JobClass, or Chain.'
                )

        elif len(args) == 2:
            # Two object arguments
            arg1, arg2 = args[0], args[1]
            result = self.data.copy()

            # Determine order and filter accordingly
            if self._is_station_or_node(arg1) and self._is_class_or_chain(arg2):
                result = self._filter_by_first_dimension_internal(result, arg1)
                result = self._filter_by_second_dimension_internal(result, arg2)
            elif self._is_class_or_chain(arg1) and self._is_station_or_node(arg2):
                result = self._filter_by_first_dimension_internal(result, arg2)
                result = self._filter_by_second_dimension_internal(result, arg1)
            else:
                raise ValueError(
                    'Expected one Station/Node and one JobClass/Chain object.'
                )

            return result

        else:
            raise ValueError('filterBy expects 1 or 2 arguments')

    def get(self, *args):
        """
        Alias for filterBy - convenient shorthand syntax.

        Args:
            *args: Same as filterBy

        Returns:
            pd.DataFrame: Filtered DataFrame
        """
        return self.filterBy(*args)

    def tget(self, *args):
        """
        Alias for filterBy - backward compatibility with tget() function.

        Args:
            *args: Same as filterBy

        Returns:
            pd.DataFrame: Filtered DataFrame
        """
        return self.filterBy(*args)

    def _is_object_arg(self, arg):
        """Check if argument is a custom object (not primitive type)."""
        # Check if it's a basic Python type
        if isinstance(arg, (int, float, str, bool, type(None))):
            return False
        if isinstance(arg, (list, tuple, dict, slice)):
            return False
        # pandas/numpy indexers are NOT model objects: a boolean mask, an Index or
        # an array must reach self.data and keep DataFrame semantics. Without this
        # they fall through to filterBy, where a pd.Series matches the
        # `hasattr(obj, 'name')` fallback -- a Series carries a .name attribute --
        # and comes back as a SILENTLY EMPTY frame. That is what made
        # `table[table['Station'] == 'Q']` return nothing.
        if isinstance(arg, (pd.Series, pd.Index, pd.DataFrame, np.ndarray)):
            return False
        # Assume everything else is an object
        return True

    def _is_station_or_node(self, obj):
        """Check if object is a Station or Node."""
        class_name = obj.__class__.__name__
        # Check explicit class names first
        if class_name in ('Station', 'Node', 'Queue', 'Sink', 'Source', 'Delay', 'Router'):
            return True
        # For duck typing support (mocks, etc.), check if object name matches first dimension columns
        if hasattr(obj, 'getName') or hasattr(obj, 'name'):
            obj_name = self._get_object_name(obj)
            columns = self.data.columns.tolist()
            if 'Station' in columns and obj_name in self.data['Station'].tolist():
                return True
            if 'Node' in columns and obj_name in self.data['Node'].tolist():
                return True
        return False

    def _is_class_or_chain(self, obj):
        """Check if object is a JobClass or Chain."""
        class_name = obj.__class__.__name__
        # Check explicit class names first
        if class_name in ('JobClass', 'OpenClass', 'ClosedClass', 'Chain'):
            return True
        # For duck typing support (mocks, etc.), check if object name matches second dimension columns
        if hasattr(obj, 'getName') or hasattr(obj, 'name'):
            obj_name = self._get_object_name(obj)
            columns = self.data.columns.tolist()
            if 'JobClass' in columns and obj_name in self.data['JobClass'].tolist():
                return True
            if 'Chain' in columns and obj_name in self.data['Chain'].tolist():
                return True
        return False

    def _filter_by_first_dimension(self, obj):
        """Filter by first dimension object (Station or Node)."""
        return self._filter_by_first_dimension_internal(self.data, obj)

    def _filter_by_first_dimension_internal(self, df, obj):
        """Internal method for filtering by first dimension."""
        obj_name = self._get_object_name(obj)
        columns = df.columns.tolist()

        # Try Station first, then Node
        if 'Station' in columns:
            return df[df['Station'] == obj_name].copy()
        elif 'Node' in columns:
            return df[df['Node'] == obj_name].copy()
        else:
            return df.copy()

    def _filter_by_second_dimension(self, obj):
        """Filter by second dimension object (JobClass or Chain)."""
        return self._filter_by_second_dimension_internal(self.data, obj)

    def _filter_by_second_dimension_internal(self, df, obj):
        """Internal method for filtering by second dimension."""
        obj_name = self._get_object_name(obj)
        columns = df.columns.tolist()

        # Try JobClass first, then Chain
        if 'JobClass' in columns:
            return df[df['JobClass'] == obj_name].copy()
        elif 'Chain' in columns:
            return df[df['Chain'] == obj_name].copy()
        else:
            return df.copy()

    def _get_object_name(self, obj):
        """Extract name from object, handling both Python and Java objects."""
        # Try Python method
        if hasattr(obj, 'getName'):
            name = obj.getName()
            return str(name) if name is not None else ''
        # Try property
        if hasattr(obj, 'name'):
            name = obj.name
            return str(name) if name is not None else ''
        # Fallback to string representation
        return str(obj)

    # Pandas compatibility methods

    def __len__(self):
        """Return number of rows."""
        return len(self.data)

    def __repr__(self):
        """Return string representation with MATLAB-style number formatting."""
        return self._format_table()

    def __str__(self):
        """Return string representation with MATLAB-style number formatting."""
        return self._format_table()

    def tabulate(self):
        """The underlying DataFrame.

        Callers written before this wrapper existed test `isinstance(t,
        DataFrame)` and fall back to `t.tabulate()`, so a table that is neither
        raises inside pandas' own `__getattr__` with a message naming DataFrame
        rather than this class. Kept as the escape hatch that contract expects;
        `.data` is the same object under its current name.
        """
        return self.data

    def to_string(self, index=False, **kwargs):
        """MATLAB-style rendering, so `print(t)` and `print(t.to_string())` agree.

        WITHOUT THIS, `to_string` fell through `__getattr__` to pandas, which
        formats to `display.precision` DECIMAL PLACES (5, set in __init__.py)
        while every other codebase prints 5 SIGNIFICANT DIGITS. The two coincide
        only for values of order 1: a cache hit rate of 0.024273 printed as
        0.02427, one digit short, and the parity comparator read a 1.2e-4
        relative gap against MATLAB's 0.024273 where the underlying values agreed
        to twelve digits. It looked like a solver defect and was a formatter.

        A caller that explicitly wants the index, or passes any other pandas
        option, gets pandas' own rendering: this override exists to fix the
        DEFAULT, not to reimplement `DataFrame.to_string`.
        """
        if index or kwargs:
            return self.data.to_string(index=index, **kwargs)
        return self._format_table()

    def _format_table(self):
        """Format the table with MATLAB-style number formatting, fitting rows
        on single lines. Name columns are left-aligned and numeric columns
        right-aligned with a two-space gap and no row index, matching the
        MATLAB and JAR table displays."""
        numeric_cols = set(self.data.select_dtypes(include=[np.number]).columns)

        columns = list(self.data.columns)
        cells = {}
        widths = {}
        for col in columns:
            if col in numeric_cols:
                vals = [_matlab_format_float(v) for v in self.data[col]]
            else:
                vals = [str(v) for v in self.data[col]]
            cells[col] = vals
            widths[col] = max([len(str(col))] + [len(v) for v in vals])

        def _row(values):
            parts = []
            for col, v in zip(columns, values):
                if col in numeric_cols:
                    parts.append(v.rjust(widths[col]))
                else:
                    parts.append(v.ljust(widths[col]))
            return '  '.join(parts).rstrip()

        lines = [_row([str(c) for c in columns])]
        for i in range(len(self.data)):
            lines.append(_row([cells[col][i] for col in columns]))
        return '\n'.join(lines)

    @property
    def shape(self):
        """Return shape of underlying DataFrame."""
        return self.data.shape

    @property
    def columns(self):
        """Return columns of underlying DataFrame."""
        return self.data.columns

    @property
    def index(self):
        """Return index of underlying DataFrame."""
        return self.data.index

    def head(self, n=5):
        """Return first n rows."""
        return self.data.head(n)

    def tail(self, n=5):
        """Return last n rows."""
        return self.data.tail(n)

    def info(self):
        """Print DataFrame info."""
        return self.data.info()

    def describe(self):
        """Return statistical description."""
        return self.data.describe()

    def to_csv(self, *args, **kwargs):
        """Export to CSV."""
        return self.data.to_csv(*args, **kwargs)

    def to_excel(self, *args, **kwargs):
        """Export to Excel."""
        return self.data.to_excel(*args, **kwargs)
