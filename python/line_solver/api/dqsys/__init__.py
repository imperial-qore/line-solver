"""
Native Python implementations for discrete-time (slotted) queueing systems.

Every function here observes the system on a lattice of unit slots rather than
on the continuous time axis, under the late-arrival rule with the departure
resolved before the arrival (Daduna's LA and D/A rules). The continuous-time
counterparts live in :mod:`line_solver.api.qsys`.

Key algorithms:
    Geo/Geo/1: dqsys_geogeo1
    Geo^X/Geo/1 (batch arrivals): dqsys_geoxgeo1, dqsys_geoxgeo1_moments
    State dependent Bernoulli server: dqsys_bernoulli1
"""

from .discrete import dqsys_geogeo1, dqsys_geoxgeo1, dqsys_geoxgeo1_moments
from .bernoulli import dqsys_bernoulli1

__all__ = [
    'dqsys_geogeo1',
    'dqsys_geoxgeo1',
    'dqsys_geoxgeo1_moments',
    'dqsys_bernoulli1',
]
