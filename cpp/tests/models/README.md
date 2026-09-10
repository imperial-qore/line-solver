# mp_pfqn model library (vendored)

The `.qn` models `cpp/tests/test_pfqn_nc_oracle.cpp` reads, copied verbatim
from `mp_pfqn.git/models/` (BSD-3, the same license as LINE). Only the models
the oracle actually names are here; the upstream library holds a few more.

They live in the tree because the oracle used to point at an absolute path in a
personal checkout, which no worker host could see: every row hit the
"file absent" branch and all three cases reported success having compared
nothing. Test data a case cannot run without belongs next to the case.

Format is mp_pfqn's own: station count, then per-class populations, demands and
(for the load-dependent models) an `MU` block. `line::io::read_qn` parses it.

To point the oracle at an upstream checkout instead, set
`LINE_MP_PFQN_MODELS=/path/to/mp_pfqn.git/models`.
