"""Session-wide fixtures for the native Python suite.

THE SESSION VERBOSITY IS PER-TEST STATE, NOT COLLECTION STATE.

Eleven test modules set ``GlobalConstants.setVerbose(VerboseLevel.SILENT)`` at
MODULE scope to keep their own solves quiet. pytest imports every collected
module before it runs anything, so in a whole-suite run those statements execute
first and leave the WHOLE SESSION silent -- including the tests that assert on
what LINE printed, which then see nothing and fail for a reason that has nothing
to do with them. Running those files alone, or in a narrow selection, hides it
entirely, which is why the pollution survived: it depends on what else was
collected.

It was invisible for another reason too, until 2026-09-11. ``api/io/logging.py``
declared a ``VerboseLevel`` of its own, a different Enum class from the one in
``line_solver.constants``, so ``printf`` and ``warning`` compared members of two
unrelated enums, never matched ``SILENT``, and printed through a silenced
session. With the two enums unified SILENT began to silence, and the pollution
became visible as four failures in ``test_mam_me_warning.py``.

The fixture restores the level before each test rather than after, because the
statement that sets it runs at IMPORT time and so precedes every test.
"""

import pytest

from line_solver.constants import GlobalConstants, VerboseLevel


@pytest.fixture(autouse=True)
def line_session_verbosity():
    """Start each test at the session default, whatever collection left behind.

    A test that wants another level sets it in its own body; only the value
    leaked by a module-scope statement in an unrelated file is undone here.
    """
    GlobalConstants.setVerbose(VerboseLevel.STD)
    yield
