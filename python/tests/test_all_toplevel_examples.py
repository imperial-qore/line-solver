"""Execute the notebooks that sit outside basic/, advanced/ and gettingstarted/.

`discrete/`, `inference/`, `opt/` and `solvers/`, plus any notebook directly in
`examples/`, had no notebooks at all before the suite became a derived artifact
and therefore no coverage. `json/` holds model-builder modules with nothing to
run, so it contributes no notebook and is deliberately absent from `roots`.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from notebook_suite import NotebookSuite, attach


class TestTopLevelExamples(NotebookSuite):
    __test__ = True
    # '' is the examples/ directory itself, collected without descending.
    roots = ['discrete', 'inference', 'opt', 'solvers', '']


attach(TestTopLevelExamples)


if __name__ == '__main__':
    import unittest
    unittest.main()
