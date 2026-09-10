"""The example notebooks must match what `tools/gen_notebooks.py` derives.

Without this, editing an example `.py` leaves its `.ipynb` behind silently, and
the notebook suites go on asserting the OLD model -- which is precisely how the
hand-authored notebooks drifted until not one of 200 matched its script.

Run `cd python && python3 tools/gen_notebooks.py` to fix a failure here.
"""

import os
import sys
import unittest

PYTHON_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
sys.path.insert(0, os.path.join(PYTHON_ROOT, 'tools'))

import gen_notebooks


class TestNotebooksAreGenerated(unittest.TestCase):

    def setUp(self):
        self._cwd = os.getcwd()
        os.chdir(PYTHON_ROOT)

    def tearDown(self):
        os.chdir(self._cwd)

    def test_every_example_has_a_notebook(self):
        missing = [p for p in gen_notebooks.find_examples()
                   if not os.path.exists(os.path.splitext(p)[0] + '.ipynb')]
        self.assertEqual(missing, [],
                         "example script(s) with no notebook; run "
                         "python3 tools/gen_notebooks.py:\n  " + "\n  ".join(missing))

    def test_no_notebook_without_a_script(self):
        """A notebook whose script is gone is unowned and unregenerable."""
        keep = set(gen_notebooks.find_examples())
        orphans = []
        for dirpath, dirnames, filenames in os.walk(gen_notebooks.EXAMPLES_DIR):
            dirnames.sort()
            for name in sorted(filenames):
                if not name.endswith('.ipynb') or '.ipynb_checkpoints' in dirpath:
                    continue
                nb = os.path.join(dirpath, name)
                if os.path.splitext(nb)[0] + '.py' not in keep:
                    orphans.append(nb)
        self.assertEqual(orphans, [],
                         "notebook(s) with no generating script; delete them or "
                         "restore the .py:\n  " + "\n  ".join(orphans))

    def test_notebooks_match_their_scripts(self):
        stale, failed = [], []
        for path in gen_notebooks.find_examples():
            nbpath = os.path.splitext(path)[0] + '.ipynb'
            try:
                text = gen_notebooks.render(path)
            except (gen_notebooks.GenError, SyntaxError) as exc:
                failed.append('%s: %s' % (path, exc))
                continue
            if not os.path.exists(nbpath):
                stale.append(nbpath)
            elif open(nbpath).read() != text:
                stale.append(nbpath)

        self.assertEqual(failed, [],
                         "script(s) the generator refuses to convert:\n  "
                         + "\n  ".join(failed))
        self.assertEqual(stale, [],
                         "notebook(s) out of date with their .py; run "
                         "python3 tools/gen_notebooks.py:\n  " + "\n  ".join(stale))


if __name__ == '__main__':
    unittest.main()
