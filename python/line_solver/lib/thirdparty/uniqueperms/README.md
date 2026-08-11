# uniqueperms

This is a derivative work ported from MATLAB to Python.

The original MATLAB implementation was written by John D'Errico
(woodchips@rochester.rr.com), Release 1.0, 2/25/08.

The Python port adds overflow protection (MAX_UNIQUE_PERMS threshold)
to prevent combinatorial explosion for large multisets, returning a
single sorted representative permutation when the count exceeds the limit.
