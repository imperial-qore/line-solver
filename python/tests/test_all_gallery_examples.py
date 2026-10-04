"""Execute every notebook under `python/examples/gallery/`.

`gallery/` was in no suite's collection list, which is why 15 of its notebooks
sat in the tree unparseable and the rest ran a CTMC the gallery scripts never
asked for. It is the single largest example folder, so it gets its own suite
rather than being folded into another.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from notebook_suite import NotebookSuite, attach


class TestGalleryExamples(NotebookSuite):
    __test__ = True
    roots = ['gallery']


attach(TestGalleryExamples)


if __name__ == '__main__':
    import unittest
    unittest.main()
