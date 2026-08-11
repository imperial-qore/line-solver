"""
API-alias resolution shared by the native solver and model classes.

Provides a __getattr__ fallback that maps the documented alternative
spellings of a capability onto the method actually implemented on the class:

- snake_case -> camelCase (get_tran_avg -> getTranAvg)
- accessor-prefixed -> bare (getAvgTable -> avg_table-style wrapper methods)
- snake_case bare -> get-prefixed camel (tran_avg -> getTranAvg)

Bare camelCase names are NOT expanded to get-prefixed forms: doing so would
make hasattr(obj, 'model') true on any class with getModel, which breaks
feature-probing code. The fallback fires only when normal lookup fails, so it
never interferes with copy/pickle or genuine AttributeError for typos.
"""

import re as _re

_CAMEL_RE = _re.compile(r'(?<!^)(?<=[a-z0-9])([A-Z])')

_ACCESSOR_PREFIXES = ('get', 'set', 'is', 'has')


# Compound metric tokens whose camel spelling is not a plain capitalization
_SPECIAL_TOKENS = {
    'respt': 'RespT', 'residt': 'ResidT', 'waitt': 'WaitT', 'passt': 'PassT',
    'qlen': 'QLen', 'sjrnt': 'SjrnT', 'arvr': 'ArvR',
    'pdf': 'PDF', 'pmf': 'PMF', 'lst': 'LST',
}


def snake_to_camel(name):
    """Convert a snake_case identifier to camelCase (get_tran_avg -> getTranAvg).

    Compound metric tokens keep their canonical spelling
    (avg_residt -> avgResidT, cdf_respt -> cdfRespT).
    """
    parts = name.split('_')
    return parts[0] + ''.join(
        _SPECIAL_TOKENS.get(p, p[:1].upper() + p[1:]) for p in parts[1:])


def camel_to_snake(name):
    """Convert a camelCase identifier to snake_case (getTranAvg -> get_tran_avg)."""
    return _CAMEL_RE.sub(r'_\1', name).lower()


def alias_candidates(name):
    """Alternative spellings under which a missing attribute may exist."""
    cands = []
    if '_' not in name and name[:1].isupper():
        # legacy PascalCase alias (GetAvgQLen -> getAvgQLen)
        camel = name[0].lower() + name[1:]
        cands.append(camel)
    else:
        camel = snake_to_camel(name) if '_' in name else name
    if camel != name and camel not in cands:
        cands.append(camel)
    snake = camel_to_snake(camel)
    if snake != name and snake not in cands:
        cands.append(snake)
    prefixed = False
    for pre in _ACCESSOR_PREFIXES:
        if camel.startswith(pre) and len(camel) > len(pre) and camel[len(pre)].isupper():
            bare = camel[len(pre):]
            bare = bare[0].lower() + bare[1:]
            cands.append(bare)
            bare_snake = camel_to_snake(bare)
            if bare_snake != bare:
                cands.append(bare_snake)
            prefixed = True
            break
    if not prefixed and ('_' in name or any(c.isupper() for c in camel[1:])):
        # see _kb/11-conventions-and-gotchas.md (Python long-tail low-hit gotchas) for rationale
        getcamel = 'get' + camel[0].upper() + camel[1:]
        cands.append(getcamel)
        get_snake = camel_to_snake(getcamel)
        if get_snake not in cands:
            cands.append(get_snake)
    return cands


def alias_getattr(self, name):
    """Shared __getattr__ fallback resolving alias spellings.

    Only invoked when normal lookup fails; raises AttributeError when the
    capability does not exist under any spelling.
    """
    if name.startswith('_'):
        raise AttributeError(name)
    cls = type(self)
    inst = object.__getattribute__(self, '__dict__')
    for cand in alias_candidates(name):
        if cand in inst:
            return inst[cand]
        if getattr(cls, cand, None) is not None:
            return getattr(self, cand)
    raise AttributeError(
        "'%s' object has no attribute '%s'" % (type(self).__name__, name))
