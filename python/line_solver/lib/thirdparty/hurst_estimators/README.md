# Hurst Parameter Estimators - Python Port

This is a derivative work ported from MATLAB to Python. The original MATLAB
implementation was written by Chu Chen (chen-chu@163.com), Version 1.0,
03/10/2008.

## Methods

The following estimators for the Hurst parameter of long-range dependent time
series are included:

- **RS** (Rescaled Range / R/S method)
- **absval** (Absolute Moment method)
- **aggvar** (Aggregate Variance method)
- **boxper** (Boxed/Modified Periodogram method)
- **diffvar** (Difference Variance method)
- **higuchi** (Higuchi's method)
- **peng** (Peng / Residuals of Regression method)
- **per** (Periodogram method)

The theories of the methods can be found in Murad S. Taqqu, Vadim Teverovsky
and Walter Willinger's paper "Estimators for long-range dependence: an
empirical study" and other related papers.

## Usage

```python
import numpy as np
from line_solver.lib.thirdparty.hurst_estimators import hurst_estimate, rs, aggvar

# Using the dispatcher
H = hurst_estimate(sequence, method='aggvar')

# Using individual estimators directly
H = rs(sequence, isplot=True)
H = aggvar(sequence)
```

## License

BSD-2-Clause. See LICENSE.txt for the full license text.
