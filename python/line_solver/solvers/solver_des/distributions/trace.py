"""
Trace file replay distribution samplers.

This module provides samplers that replay inter-arrival times or service times
from trace files, supporting various file formats.
"""

from typing import Optional, List, Union, TextIO
from pathlib import Path
import numpy as np

from .base import DistributionSampler


class FileTraceSampler(DistributionSampler):
    """
    Sampler that replays values from a trace file.

    Supports multiple file formats:
    - CSV (comma-separated)
    - TSV (tab-separated)
    - Text (newline-separated)
    - NumPy (.npy or .npz)

    The trace file should contain numeric values representing:
    - Inter-arrival times (for arrival processes)
    - Service times (for service processes)

    Attributes:
        values: Loaded trace values
        loop: Whether to loop when exhausted
        index: Current position in trace
    """

    def __init__(
        self,
        file_path: Union[str, Path],
        column: int = 0,
        skip_header: bool = False,
        loop: bool = True,
        delimiter: Optional[str] = None,
        rng: Optional[np.random.Generator] = None
    ):
        """
        Initialize from trace file.

        Args:
            file_path: Path to trace file
            column: Column index for multi-column files (default: 0)
            skip_header: Skip first line if header (default: False)
            loop: Loop back to start when exhausted (default: True)
            delimiter: Column delimiter (auto-detected if None)
            rng: Random number generator (not used for replay)
        """
        super().__init__(rng)

        self.file_path = Path(file_path)
        self.loop = loop
        self.index = 0

        # Load values from file
        self.values = self._load_trace(column, skip_header, delimiter)

        if len(self.values) == 0:
            raise ValueError(f"No valid values found in trace file: {file_path}")

    def _load_trace(
        self,
        column: int,
        skip_header: bool,
        delimiter: Optional[str]
    ) -> np.ndarray:
        """
        Load trace values from file.

        Args:
            column: Column index
            skip_header: Whether to skip header line
            delimiter: Column delimiter

        Returns:
            Array of trace values
        """
        suffix = self.file_path.suffix.lower()

        # NumPy binary format
        if suffix == '.npy':
            return np.load(self.file_path).flatten()

        if suffix == '.npz':
            data = np.load(self.file_path)
            # Use first array in the archive
            key = list(data.keys())[0]
            return data[key].flatten()

        # Text formats
        if delimiter is None:
            if suffix == '.csv':
                delimiter = ','
            elif suffix in ('.tsv', '.tab'):
                delimiter = '\t'
            else:
                # Auto-detect
                delimiter = self._detect_delimiter()

        values = []
        with open(self.file_path, 'r') as f:
            if skip_header:
                next(f)

            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                if delimiter:
                    parts = line.split(delimiter)
                    if column < len(parts):
                        try:
                            values.append(float(parts[column].strip()))
                        except ValueError:
                            continue
                else:
                    try:
                        values.append(float(line))
                    except ValueError:
                        continue

        return np.array(values)

    def _detect_delimiter(self) -> Optional[str]:
        """Detect delimiter from first line of file."""
        with open(self.file_path, 'r') as f:
            line = f.readline()

        if ',' in line:
            return ','
        elif '\t' in line:
            return '\t'
        elif ';' in line:
            return ';'
        elif ' ' in line:
            return None  # Split on whitespace
        return None

    def sample(self) -> float:
        """Get next value from trace."""
        if self.index >= len(self.values):
            if self.loop:
                self.index = 0
            else:
                return self.values[-1]

        value = self.values[self.index]
        self.index += 1
        return float(value)

    def sample_n(self, n: int) -> np.ndarray:
        """Get next n values from trace."""
        result = np.zeros(n)
        for i in range(n):
            result[i] = self.sample()
        return result

    @property
    def mean(self) -> float:
        """Mean of trace values."""
        return float(np.mean(self.values))

    @property
    def variance(self) -> float:
        """Variance of trace values."""
        return float(np.var(self.values))

    def reset(self) -> None:
        """Reset to beginning of trace."""
        self.index = 0

    @classmethod
    def from_csv(
        cls,
        file_path: Union[str, Path],
        column: int = 0,
        skip_header: bool = True,
        loop: bool = True,
        rng: Optional[np.random.Generator] = None
    ) -> 'FileTraceSampler':
        """
        Create sampler from CSV file.

        Args:
            file_path: Path to CSV file
            column: Column index (default: 0)
            skip_header: Skip header row (default: True)
            loop: Loop when exhausted (default: True)
            rng: Random number generator

        Returns:
            FileTraceSampler instance
        """
        return cls(
            file_path=file_path,
            column=column,
            skip_header=skip_header,
            loop=loop,
            delimiter=',',
            rng=rng
        )

    @classmethod
    def from_list(
        cls,
        values: List[float],
        loop: bool = True,
        rng: Optional[np.random.Generator] = None
    ) -> 'FileTraceSampler':
        """
        Create sampler from list of values (in-memory trace).

        Args:
            values: List of values
            loop: Loop when exhausted (default: True)
            rng: Random number generator

        Returns:
            FileTraceSampler instance
        """
        sampler = object.__new__(cls)
        sampler.rng = rng if rng is not None else np.random.default_rng()
        sampler.values = np.array(values)
        sampler.loop = loop
        sampler.index = 0
        sampler.file_path = None
        return sampler


class MultiColumnTraceSampler:
    """
    Sampler that replays multiple correlated columns from a trace file.

    Useful when inter-arrival times and service times are correlated
    (e.g., job size correlates with processing time).

    Each row in the file represents one job, with columns for:
    - Inter-arrival time
    - Service time(s) at different stations
    - Job class (optional)
    """

    def __init__(
        self,
        file_path: Union[str, Path],
        iat_column: int = 0,
        service_columns: List[int] = None,
        class_column: Optional[int] = None,
        skip_header: bool = True,
        delimiter: str = ',',
        loop: bool = True
    ):
        """
        Initialize from multi-column trace file.

        Args:
            file_path: Path to trace file
            iat_column: Column for inter-arrival times
            service_columns: Columns for service times at each station
            class_column: Column for job class (optional)
            skip_header: Skip header row
            delimiter: Column delimiter
            loop: Loop when exhausted
        """
        self.file_path = Path(file_path)
        self.loop = loop
        self.index = 0

        # Load all data
        self.data = self._load_data(skip_header, delimiter)

        # Extract columns
        self.iat_column = iat_column
        self.service_columns = service_columns or [1]
        self.class_column = class_column

        self.n_rows = len(self.data)
        if self.n_rows == 0:
            raise ValueError(f"No valid data found in trace file: {file_path}")

    def _load_data(
        self,
        skip_header: bool,
        delimiter: str
    ) -> np.ndarray:
        """Load all columns from file."""
        rows = []
        with open(self.file_path, 'r') as f:
            if skip_header:
                next(f)

            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                try:
                    values = [float(x.strip()) for x in line.split(delimiter)]
                    rows.append(values)
                except ValueError:
                    continue

        if not rows:
            return np.array([])

        # Pad rows to same length
        max_cols = max(len(row) for row in rows)
        for row in rows:
            while len(row) < max_cols:
                row.append(0.0)

        return np.array(rows)

    def get_next(self) -> dict:
        """
        Get next trace entry.

        Returns:
            Dictionary with:
                iat: Inter-arrival time
                service_times: List of service times
                job_class: Job class (if available)
        """
        if self.index >= self.n_rows:
            if self.loop:
                self.index = 0
            else:
                self.index = self.n_rows - 1

        row = self.data[self.index]
        self.index += 1

        result = {
            'iat': row[self.iat_column],
            'service_times': [row[c] for c in self.service_columns],
        }

        if self.class_column is not None:
            result['job_class'] = int(row[self.class_column])

        return result

    def get_iat(self) -> float:
        """Get next inter-arrival time (advances index)."""
        entry = self.get_next()
        return entry['iat']

    def get_service_time(self, station_idx: int = 0) -> float:
        """
        Get service time for a station.

        Note: Does NOT advance index. Call get_next() first.
        """
        if self.index == 0:
            row = self.data[0]
        else:
            row = self.data[self.index - 1]

        if station_idx < len(self.service_columns):
            col = self.service_columns[station_idx]
            return row[col]
        return 0.0

    def reset(self) -> None:
        """Reset to beginning of trace."""
        self.index = 0


class OnlineTraceSampler(DistributionSampler):
    """
    Sampler that reads values from a stream/pipe in real-time.

    Useful for:
    - Real-time trace replay
    - External workload generators
    - Named pipes (FIFOs)
    """

    def __init__(
        self,
        stream: TextIO,
        default_value: float = 1.0,
        rng: Optional[np.random.Generator] = None
    ):
        """
        Initialize from stream.

        Args:
            stream: Input stream (file-like object)
            default_value: Value to return if stream is empty
            rng: Random number generator
        """
        super().__init__(rng)
        self.stream = stream
        self.default_value = default_value
        self.buffer: List[float] = []

    def sample(self) -> float:
        """Get next value from stream."""
        # Try to read from buffer first
        if self.buffer:
            return self.buffer.pop(0)

        # Read from stream
        try:
            line = self.stream.readline()
            if line:
                return float(line.strip())
        except (ValueError, IOError):
            pass

        return self.default_value

    def buffer_values(self, n: int) -> None:
        """Pre-buffer n values from stream."""
        for _ in range(n):
            try:
                line = self.stream.readline()
                if line:
                    self.buffer.append(float(line.strip()))
            except (ValueError, IOError):
                break

    @property
    def mean(self) -> float:
        """Not available for online stream."""
        return self.default_value

    @property
    def variance(self) -> float:
        """Not available for online stream."""
        return 0.0


def load_trace(
    file_path: Union[str, Path],
    **kwargs
) -> FileTraceSampler:
    """
    Convenience function to load a trace file.

    Auto-detects file format and creates appropriate sampler.

    Args:
        file_path: Path to trace file
        **kwargs: Additional arguments for FileTraceSampler

    Returns:
        FileTraceSampler instance
    """
    return FileTraceSampler(file_path, **kwargs)
