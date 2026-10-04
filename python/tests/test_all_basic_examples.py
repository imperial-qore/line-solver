"""Execute every notebook under `python/examples/basic/`.

The harness lives in `notebook_suite.py`; this file only names the tree it owns.
`roots` is collected recursively, so a new `basic/<area>` folder is picked up
without editing anything here.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from notebook_suite import NotebookSuite, attach


class TestBasicExamples(NotebookSuite):
    __test__ = True
    roots = ['basic']
    # The basic models are cheap; a hung kernel is a defect, so let it hang
    # visibly rather than be cut off at an arbitrary deadline.
    timeout = None


attach(TestBasicExamples)


if __name__ == '__main__':
    import unittest
    unittest.main()
