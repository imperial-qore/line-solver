"""Shared helpers for state-dependent round-robin routing pointers.

Weighted round-robin (WRROBIN) stores its per-node pointer as a POSITION in a
weighted-outlink cycle rather than as a destination node index, because a
destination may repeat within one cycle (e.g. weights 2:1 give the cycle
[d1, d1, d2]). A node-index pointer cannot disambiguate the two occurrences of
d1, so the position representation is required. This mirrors MATLAB
refreshRoutingMatrix, which enumerates positions and maps them back to
destinations via sub_wrr.
"""


def wrr_weighted_outlinks(sn, ind, r):
    """Return the weighted round-robin cycle for node ``ind``, class ``r``.

    The cycle is a list of destination node indices with each destination
    repeated by its integer weight, as populated by
    Network._refresh_statedep_routing_params. Returns ``None`` when the node
    has no WRROBIN cycle recorded.
    """
    nodeparam = getattr(sn, 'nodeparam', None)
    if nodeparam is None or ind not in nodeparam:
        return None
    entry = nodeparam[ind]
    if not isinstance(entry, dict):
        return None
    class_entry = entry.get(r)
    if class_entry is None:
        class_entry = entry.get(int(r))
    if isinstance(class_entry, dict):
        wol = class_entry.get('weighted_outlinks')
        if wol is not None and len(wol) > 0:
            return list(wol)
    return None
