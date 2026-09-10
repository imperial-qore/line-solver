"""Execute every notebook under `python/examples/advanced/`.

The harness lives in `notebook_suite.py`; this file only names the tree it owns.
`roots` is collected recursively, so a new `advanced/<area>` folder is picked up
without editing anything here.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from notebook_suite import NotebookSuite, attach


class TestAdvancedExamples(NotebookSuite):
    __test__ = True
    roots = ['advanced']


attach(TestAdvancedExamples)


if __name__ == '__main__':
    import unittest
    unittest.main()
