"""
Factory function for creating scheduling strategy instances.

This module provides the create_scheduler function that instantiates
the appropriate scheduler based on the scheduling strategy enum.
"""

from typing import Optional, List, Any
from enum import IntEnum

from .base import SchedulingStrategy
from .fcfs import FCFSScheduler, PriorityFCFSScheduler, EDDScheduler
from .lcfs import LCFSScheduler, LCFSPriorityScheduler
from .siro import SIROScheduler
from .job_based import SJFScheduler, LJFScheduler, SEPTScheduler, LEPTScheduler


class SchedStrategyID(IntEnum):
    """Scheduling strategy identifiers."""
    FCFS = 0
    LCFS = 1
    LCFSPR = 2
    PS = 3
    DPS = 4
    GPS = 5
    INF = 6
    RAND = 7
    HOL = 8
    SEPT = 9
    LEPT = 10
    SIRO = 11
    SJF = 12
    LJF = 13
    POLLING = 14
    EXT = 15
    LCFSPI = 16
    SRPT = 17
    SRPTPRIO = 18
    EDD = 19
    EDF = 20
    PSPRIO = 21
    DPSPRIO = 22
    GPSPRIO = 23
    FCFSPRIO = 24
    LCFSPRIO = 25
    LCFSPRPRIO = 26
    LCFSPIPRIO = 27
    FCFSPR = 28
    FCFSPI = 29
    FCFSPRPRIO = 30
    FCFSPIPRIO = 31
    FB = 32
    FBPRIO = 33
    SETF = 34
    SETFPRIO = 35


def get_sched_strategy_id(strategy: Any) -> int:
    """
    Get numeric ID for a scheduling strategy.

    Args:
        strategy: SchedStrategy enum or string

    Returns:
        Integer ID for the strategy
    """
    if isinstance(strategy, int):
        return strategy

    # Handle JPype enum
    if hasattr(strategy, 'name'):
        name = str(strategy.name()).upper()
    elif hasattr(strategy, 'value'):
        name = str(strategy.value).upper()
    else:
        name = str(strategy).upper()

    # Map names to IDs
    name_to_id = {
        'FCFS': SchedStrategyID.FCFS,
        'LCFS': SchedStrategyID.LCFS,
        'LCFSPR': SchedStrategyID.LCFSPR,
        'PS': SchedStrategyID.PS,
        'DPS': SchedStrategyID.DPS,
        'GPS': SchedStrategyID.GPS,
        'INF': SchedStrategyID.INF,
        'RAND': SchedStrategyID.RAND,
        'HOL': SchedStrategyID.HOL,
        'SEPT': SchedStrategyID.SEPT,
        'LEPT': SchedStrategyID.LEPT,
        'SIRO': SchedStrategyID.SIRO,
        'SJF': SchedStrategyID.SJF,
        'LJF': SchedStrategyID.LJF,
        'POLLING': SchedStrategyID.POLLING,
        'EXT': SchedStrategyID.EXT,
        'LCFSPI': SchedStrategyID.LCFSPI,
        'SRPT': SchedStrategyID.SRPT,
        'SRPTPRIO': SchedStrategyID.SRPTPRIO,
        'EDD': SchedStrategyID.EDD,
        'EDF': SchedStrategyID.EDF,
        'PSPRIO': SchedStrategyID.PSPRIO,
        'DPSPRIO': SchedStrategyID.DPSPRIO,
        'GPSPRIO': SchedStrategyID.GPSPRIO,
        'FCFSPRIO': SchedStrategyID.FCFSPRIO,
        'LCFSPRIO': SchedStrategyID.LCFSPRIO,
        'LCFSPRPRIO': SchedStrategyID.LCFSPRPRIO,
        'LCFSPIPRIO': SchedStrategyID.LCFSPIPRIO,
        'FCFSPR': SchedStrategyID.FCFSPR,
        'FCFSPI': SchedStrategyID.FCFSPI,
        'FCFSPRPRIO': SchedStrategyID.FCFSPRPRIO,
        'FCFSPIPRIO': SchedStrategyID.FCFSPIPRIO,
        'FB': SchedStrategyID.FB,
        'LAS': SchedStrategyID.FB,  # Alias: Least Attained Service
        'FBPRIO': SchedStrategyID.FBPRIO,
        'SETF': SchedStrategyID.SETF,
        'SETFPRIO': SchedStrategyID.SETFPRIO,
    }

    return name_to_id.get(name, SchedStrategyID.FCFS)


def create_scheduler(
    strategy: Any,
    num_classes: int,
    num_servers: int,
    service_rates: Optional[List[float]] = None,
    class_weights: Optional[List[float]] = None,
    class_priorities: Optional[List[int]] = None,
    **kwargs
) -> SchedulingStrategy:
    """
    Create a scheduling strategy instance.

    Args:
        strategy: SchedStrategy enum or ID
        num_classes: Number of job classes
        num_servers: Number of servers at the station
        service_rates: Service rates per class (for SEPT/LEPT)
        class_weights: Weights per class (for DPS/GPS)
        class_priorities: Priorities per class

    Returns:
        SchedulingStrategy instance

    Raises:
        ValueError: If strategy is not supported
    """
    strategy_id = get_sched_strategy_id(strategy)

    # Non-preemptive strategies
    if strategy_id == SchedStrategyID.FCFS:
        return FCFSScheduler(num_classes, num_servers)

    elif strategy_id in (SchedStrategyID.HOL, SchedStrategyID.FCFSPRIO):
        return PriorityFCFSScheduler(num_classes, num_servers)

    elif strategy_id == SchedStrategyID.LCFS:
        return LCFSScheduler(num_classes, num_servers)

    elif strategy_id == SchedStrategyID.LCFSPRIO:
        return LCFSPriorityScheduler(num_classes, num_servers)

    elif strategy_id == SchedStrategyID.SIRO:
        return SIROScheduler(num_classes, num_servers)

    elif strategy_id == SchedStrategyID.SJF:
        return SJFScheduler(num_classes, num_servers)

    elif strategy_id == SchedStrategyID.LJF:
        return LJFScheduler(num_classes, num_servers)

    elif strategy_id == SchedStrategyID.SEPT:
        return SEPTScheduler(num_classes, num_servers, service_rates)

    elif strategy_id == SchedStrategyID.LEPT:
        return LEPTScheduler(num_classes, num_servers, service_rates)

    elif strategy_id == SchedStrategyID.EDD:
        return EDDScheduler(num_classes, num_servers)

    # Preemptive strategies - import here to avoid circular imports
    elif strategy_id in (SchedStrategyID.LCFSPR, SchedStrategyID.LCFSPRPRIO):
        from .preemptive import LCFSPRScheduler
        has_priority = (strategy_id == SchedStrategyID.LCFSPRPRIO)
        return LCFSPRScheduler(num_classes, num_servers, has_priority)

    elif strategy_id in (SchedStrategyID.LCFSPI, SchedStrategyID.LCFSPIPRIO):
        from .preemptive import LCFSPIScheduler
        has_priority = (strategy_id == SchedStrategyID.LCFSPIPRIO)
        return LCFSPIScheduler(num_classes, num_servers, has_priority)

    elif strategy_id in (SchedStrategyID.FCFSPR, SchedStrategyID.FCFSPRPRIO):
        from .preemptive import FCFSPRScheduler
        has_priority = (strategy_id == SchedStrategyID.FCFSPRPRIO)
        return FCFSPRScheduler(num_classes, num_servers, has_priority)

    elif strategy_id in (SchedStrategyID.FCFSPI, SchedStrategyID.FCFSPIPRIO):
        from .preemptive import FCFSPIScheduler
        has_priority = (strategy_id == SchedStrategyID.FCFSPIPRIO)
        return FCFSPIScheduler(num_classes, num_servers, has_priority)

    elif strategy_id in (SchedStrategyID.SRPT, SchedStrategyID.SRPTPRIO):
        from .preemptive import SRPTScheduler
        has_priority = (strategy_id == SchedStrategyID.SRPTPRIO)
        return SRPTScheduler(num_classes, num_servers, has_priority)

    elif strategy_id == SchedStrategyID.EDF:
        from .preemptive import EDFScheduler
        return EDFScheduler(num_classes, num_servers)

    elif strategy_id in (SchedStrategyID.FB, SchedStrategyID.FBPRIO):
        from .preemptive import FBScheduler
        has_priority = (strategy_id == SchedStrategyID.FBPRIO)
        return FBScheduler(num_classes, num_servers, has_priority)

    elif strategy_id in (SchedStrategyID.SETF, SchedStrategyID.SETFPRIO):
        from .preemptive import SETFScheduler
        has_priority = (strategy_id == SchedStrategyID.SETFPRIO)
        return SETFScheduler(num_classes, num_servers, has_priority)

    # Processor sharing strategies
    elif strategy_id in (SchedStrategyID.PS, SchedStrategyID.PSPRIO):
        from .ps import PSScheduler
        has_priority = (strategy_id == SchedStrategyID.PSPRIO)
        return PSScheduler(num_classes, num_servers, has_priority)

    elif strategy_id in (SchedStrategyID.DPS, SchedStrategyID.DPSPRIO):
        from .ps import DPSScheduler
        has_priority = (strategy_id == SchedStrategyID.DPSPRIO)
        return DPSScheduler(num_classes, num_servers, class_weights, has_priority)

    elif strategy_id in (SchedStrategyID.GPS, SchedStrategyID.GPSPRIO):
        from .ps import GPSScheduler
        has_priority = (strategy_id == SchedStrategyID.GPSPRIO)
        return GPSScheduler(num_classes, num_servers, class_weights, has_priority)

    # Polling
    elif strategy_id == SchedStrategyID.POLLING:
        from .polling import PollingScheduler
        return PollingScheduler(num_classes, num_servers, **kwargs)

    # Infinite server (delay node) - use FCFS with infinite servers
    elif strategy_id == SchedStrategyID.INF:
        return FCFSScheduler(num_classes, num_servers)

    # Default to FCFS
    else:
        return FCFSScheduler(num_classes, num_servers)
