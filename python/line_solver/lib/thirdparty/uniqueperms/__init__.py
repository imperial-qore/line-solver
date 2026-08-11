"""
uniqueperms - Generate all unique permutations of a vector with duplicate elements.

Ported from MATLAB to Python. Original MATLAB implementation by:
    Author: John D'Errico
    e-mail: woodchips@rochester.rr.com
    Release: 1.0
    Release date: 2/25/08
"""

from .uniqueperms import uniqueperms, uniqueperms_recursive, MAX_UNIQUE_PERMS

__all__ = ["uniqueperms", "uniqueperms_recursive", "MAX_UNIQUE_PERMS"]
