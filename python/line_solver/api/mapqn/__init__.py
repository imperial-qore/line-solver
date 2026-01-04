"""
MAP Queueing Network (MAPQN) Bounds Computation.

Native Python implementations for computing performance bounds in MAP
queueing networks using linear programming formulations.

Key algorithms:
    mapqn_bnd_lr_pf: Linear reduction bounds for product-form networks
    mapqn_bnd_lr: General linear reduction bounds with phases
    mapqn_bnd_lr_mva: MVA-based linear reduction bounds
    mapqn_bnd_qr: General quadratic reduction bounds
    mapqn_bnd_qr_delay: Quadratic reduction bounds for delay systems
    mapqn_bnd_qr_ld: Quadratic reduction bounds for load-dependent systems
    mapqn_qr_bounds_bas: QR bounds with blocking-after-service
    mapqn_qr_bounds_rsrd: QR bounds with RSRD blocking

Key classes:
    MapqnSolution: Container for LP optimization results
    MapqnLpModel: Base LP model builder
    PFParameters: Product-form model parameters
    LinearReductionParameters: General linear reduction parameters

References:
    Casale, G., et al. "LINE: A unified library for queueing network modeling."
"""

from .solution import MapqnSolution

from .lpmodel import (
    MapqnLpModel,
    LinearConstraint,
    LinearConstraintBuilder,
)

from .parameters import (
    MapqnParameters,
    PFParameters,
    LinearReductionParameters,
    MVAVersionParameters,
    QRBoundsBasParameters,
    QRBoundsRsrdParameters,
    QuadraticDelayParameters,
    QuadraticLDParameters,
)

from .bnd_lr_pf import mapqn_bnd_lr_pf
from .bnd_lr import mapqn_bnd_lr
from .bnd_lr_mva import mapqn_bnd_lr_mva
from .bnd_qr import mapqn_bnd_qr
from .bnd_qr_delay import mapqn_bnd_qr_delay
from .bnd_qr_ld import mapqn_bnd_qr_ld
from .qr_bounds_bas import mapqn_qr_bounds_bas
from .qr_bounds_rsrd import mapqn_qr_bounds_rsrd


__all__ = [
    # Solution container
    'MapqnSolution',
    # LP model infrastructure
    'MapqnLpModel',
    'LinearConstraint',
    'LinearConstraintBuilder',
    # Parameter classes
    'MapqnParameters',
    'PFParameters',
    'LinearReductionParameters',
    'MVAVersionParameters',
    'QRBoundsBasParameters',
    'QRBoundsRsrdParameters',
    'QuadraticDelayParameters',
    'QuadraticLDParameters',
    # Bound algorithms - Linear Reduction
    'mapqn_bnd_lr_pf',
    'mapqn_bnd_lr',
    'mapqn_bnd_lr_mva',
    # Bound algorithms - Quadratic Reduction
    'mapqn_bnd_qr',
    'mapqn_bnd_qr_delay',
    'mapqn_bnd_qr_ld',
    # Bound algorithms - Blocking
    'mapqn_qr_bounds_bas',
    'mapqn_qr_bounds_rsrd',
]
