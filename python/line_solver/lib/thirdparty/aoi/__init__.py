"""
Age of Information (AoI) Analysis.

This module implements comprehensive AoI analysis tools including:
- Analytical formulas for standard queues (M/M/1, M/G/1, G/M/1, etc.)
- Laplace-Stieltjes Transform (LST) functions for distributions
- Markovian fluid queue (MFQ) solvers for general systems
- Network topology validation and parameter extraction

**Scheduling Disciplines**:
- FCFS: First-Come First-Served
- LCFS-PR: Last-Come First-Served with Preemption
- LCFS-S: LCFS with Service discarding
- LCFS-D: LCFS with Departure discarding

**Supported Topologies** (MFQ solvers):
- Bufferless (capacity=1): PH/PH/1/1 or PH/PH/1/1* (preemptive service)
- Single-buffer (capacity=2): M/PH/1/2 or M/PH/1/2* (replacement service)

**References**:
- Y. Inoue, H. Masuyama, T. Takine, T. Tanaka, "A General Formula for
  the Stationary Distribution of the Age of Information and Its
  Application to Single-Server Queues," IEEE Trans. IT, 2019.
- Dogan, O., Akar, N., & Atay, F. F., "Age of Information in Markovian
  Fluid Queues", arXiv:2003.09408, 2020.

**Key Functions**:
- LST functions: aoi_lst_exp, aoi_lst_erlang, aoi_lst_det, aoi_lst_ph
- FCFS queues: aoi_fcfs_mm1, aoi_fcfs_md1, aoi_fcfs_dm1, aoi_fcfs_mgi1, aoi_fcfs_gim1
- LCFS queues: aoi_lcfspr_mm1, aoi_lcfspr_mgi1, aoi_lcfspr_gim1, etc.
- MFQ solvers: solve_bufferless, solve_singlebuffer
- Validation: aoi_is_aoi, aoi_extract_params
"""

from .validation import aoi_is_aoi, aoi_extract_params
from .conversion import aoi_dist2ph
from .bufferless import solve_bufferless
from .singlebuffer import solve_singlebuffer

from .lst import (
    aoi_lst_exp,
    aoi_lst_erlang,
    aoi_lst_det,
    aoi_lst_ph,
)

from .analytical import (
    # FCFS queues
    aoi_fcfs_mm1,
    aoi_fcfs_md1,
    aoi_fcfs_dm1,
    aoi_fcfs_mgi1,
    aoi_fcfs_gim1,
    # LCFS preemptive queues
    aoi_lcfspr_mm1,
    aoi_lcfspr_md1,
    aoi_lcfspr_dm1,
    aoi_lcfspr_mgi1,
    aoi_lcfspr_gim1,
    # LCFS with discarding
    aoi_lcfss_mgi1,
    aoi_lcfss_gim1,
    aoi_lcfsd_mgi1,
    aoi_lcfsd_gim1,
)

__all__ = [
    # Validation and extraction
    'aoi_is_aoi',
    'aoi_extract_params',
    'aoi_dist2ph',
    # MFQ solvers
    'solve_bufferless',
    'solve_singlebuffer',
    # LST functions
    'aoi_lst_exp',
    'aoi_lst_erlang',
    'aoi_lst_det',
    'aoi_lst_ph',
    # FCFS queues
    'aoi_fcfs_mm1',
    'aoi_fcfs_md1',
    'aoi_fcfs_dm1',
    'aoi_fcfs_mgi1',
    'aoi_fcfs_gim1',
    # LCFS preemptive queues
    'aoi_lcfspr_mm1',
    'aoi_lcfspr_md1',
    'aoi_lcfspr_dm1',
    'aoi_lcfspr_mgi1',
    'aoi_lcfspr_gim1',
    # LCFS with discarding
    'aoi_lcfss_mgi1',
    'aoi_lcfss_gim1',
    'aoi_lcfsd_mgi1',
    'aoi_lcfsd_gim1',
]
