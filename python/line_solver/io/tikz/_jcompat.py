"""
Java-compatible primitives shared by the TikZ exporter.

The TikZ port is meant to emit the same text as the JAR for the same model, so
number formatting, identifier sanitization and min/max follow Java semantics.
"""

import math
from decimal import Decimal, ROUND_HALF_UP

# java.lang.Double.MIN_VALUE is the smallest POSITIVE double, and the JAR seeds its running maxima with it
JAVA_DOUBLE_MIN_VALUE = 5e-324
JAVA_DOUBLE_MAX_VALUE = 1.7976931348623157e308


def jformat(value, digits):
    """Format a double as Java's ``String.format("%.<digits>f", value)`` does.

    Java rounds HALF_UP on the shortest decimal representation of the double,
    whereas Python's ``'%.2f'`` rounds the exact binary value half-to-even, so
    e.g. 0.125 gives "0.13" in Java and "0.12" in Python.
    """
    value = float(value)
    if math.isnan(value):
        return 'NaN'
    if math.isinf(value):
        return 'Infinity' if value > 0 else '-Infinity'
    quantum = Decimal(1).scaleb(-digits)
    return str(Decimal(repr(value)).quantize(quantum, rounding=ROUND_HALF_UP))


def f2(value):
    """Java ``%.2f``."""
    return jformat(value, 2)


def jmax(a, b):
    """java.lang.Math.max for non-NaN doubles (0.0 beats -0.0)."""
    if a == b == 0.0:
        return b if math.copysign(1.0, a) < 0 else a
    return a if a >= b else b


def jmin(a, b):
    """java.lang.Math.min for non-NaN doubles (-0.0 beats 0.0)."""
    if a == b == 0.0:
        return a if math.copysign(1.0, a) < 0 else b
    return a if a <= b else b


def sanitize_id(name):
    """Java ``name.replaceAll("[^a-zA-Z0-9]", "_")``; a non-BMP character is two UTF-16 units, so two underscores."""
    out = []
    for ch in name:
        if ('a' <= ch <= 'z') or ('A' <= ch <= 'Z') or ('0' <= ch <= '9'):
            out.append(ch)
        elif ord(ch) > 0xFFFF:
            out.append('__')
        else:
            out.append('_')
    return ''.join(out)


def escape_latex(text):
    """Escape LaTeX special characters, applying the replacements in the JAR's order."""
    return (text
            .replace('\\', '\\textbackslash{}')
            .replace('_', '\\_')
            .replace('&', '\\&')
            .replace('%', '\\%')
            .replace('$', '\\$')
            .replace('#', '\\#')
            .replace('{', '\\{')
            .replace('}', '\\}')
            .replace('~', '\\textasciitilde{}')
            .replace('^', '\\textasciicircum{}'))


def node_kind(node):
    """Classify a node the way the JAR's instanceof chain does (Delay before Queue, since Delay extends Queue)."""
    from ...lang import nodes as N
    for kind, cls in (('delay', N.Delay), ('queue', N.Queue), ('source', N.Source), ('sink', N.Sink),
                      ('fork', N.Fork), ('join', N.Join), ('router', N.Router), ('classswitch', N.ClassSwitch),
                      ('cache', N.Cache), ('logger', N.Logger), ('place', N.Place), ('transition', N.Transition)):
        if isinstance(node, cls):
            return kind
    return 'generic'
