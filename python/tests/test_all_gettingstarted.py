"""Execute every notebook under `python/examples/gettingstarted/`.

This suite used to list its tutorials by hand and stopped at tut10, so
tut11-tut13 ran nowhere and tut11's `TypeError` in `SolverLN._ph_delay_mean`
went unseen. It now discovers the folder, so a new tutorial is covered the
moment its notebook is generated.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from notebook_suite import NotebookSuite, attach


class TestGettingStarted(NotebookSuite):
    __test__ = True
    roots = ['gettingstarted']


attach(TestGettingStarted)


if __name__ == '__main__':
    import unittest
    unittest.main()
