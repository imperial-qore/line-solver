"""Union of every notebook example suite: one entry point for the whole tree.

`run-tests.sh` --ignore's this file, because collecting it alongside the five
per-area suites would run every notebook twice. It stays as the manual
whole-suite entry point (`./run-tests.sh -f test_all_examples.py`) and as the
guard that no example folder belongs to no suite: `test_every_notebook_is_owned`
fails if a notebook this file finds is missing from all the per-area suites.
"""

import os
import sys
import unittest

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from notebook_suite import NotebookSuite, attach


class TestExamples(NotebookSuite):
    __test__ = True
    # '' is the examples/ directory itself, collected without descending.
    roots = ['basic', 'advanced', 'gettingstarted', 'gallery',
             'discrete', 'inference', 'opt', 'solvers', '']

    def test_every_notebook_is_owned_by_a_suite(self):
        """No notebook may sit outside the per-area suites.

        A hardcoded collection list is how 72 notebooks -- the whole of
        gallery/, basic/workflowModels, advanced/agentModels, advanced/fcRegion
        and gettingstarted tut11-tut13 -- ended up running nowhere. This asserts
        the union of the per-area suites still covers everything.
        """
        from test_all_advanced_examples import TestAdvancedExamples
        from test_all_basic_examples import TestBasicExamples
        from test_all_gallery_examples import TestGalleryExamples
        from test_all_gettingstarted import TestGettingStarted
        from test_all_toplevel_examples import TestTopLevelExamples

        owned = set()
        for suite in (TestBasicExamples, TestAdvancedExamples, TestGettingStarted,
                      TestGalleryExamples, TestTopLevelExamples):
            owned.update(suite.discover_notebooks())

        orphans = sorted(set(self.discover_notebooks()) - owned)
        self.assertEqual(orphans, [],
                         "notebook(s) in no per-area suite; add a root to one of "
                         "the test_all_*_examples.py suites:\n  "
                         + "\n  ".join(orphans))


attach(TestExamples)


if __name__ == '__main__':
    unittest.main()
