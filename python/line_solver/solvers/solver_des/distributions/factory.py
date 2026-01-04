"""
Distribution sampler factory for creating samplers from NetworkStruct.

This module provides factory functions to create distribution samplers
based on the model specification in NetworkStruct, supporting both
rate-based and explicit distribution parameters.
"""

from typing import Dict, Optional, Tuple, Any
import numpy as np

from .base import DistributionSampler, ImmediateSampler, DisabledSampler
from .continuous import (
    ExponentialSampler,
    ErlangSampler,
    GammaSampler,
    UniformSampler,
    DeterministicSampler,
    ParetoSampler,
    WeibullSampler,
    LognormalSampler,
    NormalSampler,
    BetaSampler,
    MixtureSampler,
    ReplayTraceSampler,
)
from .phase_type import (
    PHSampler,
    APHSampler,
    CoxianSampler,
    Coxian2Sampler,
    HyperExponentialSampler,
    HyperExp2Sampler,
)
from .map_bmap import (
    MAPSampler,
    MAP2Sampler,
    BMAPSampler,
    MMPPSampler,
)
from .discrete import (
    PoissonSampler,
    BinomialSampler,
    GeometricSampler,
    NegativeBinomialSampler,
    BernoulliSampler,
    DiscreteUniformSampler,
    ZipfSampler,
    CategoricalSampler,
    EmpiricalDiscreteSampler,
)
from .trace import (
    FileTraceSampler,
    MultiColumnTraceSampler,
)


class DistributionFactory:
    """
    Factory for creating distribution samplers from NetworkStruct.

    The factory supports multiple input formats:
    1. Rate-based (scalar rate parameter) -> Exponential
    2. Explicit distribution types (EXP, ERLANG, PH, MAP, etc.)
    3. Phase-type parameters (moments, coefficients, etc.)
    4. MAP/BMAP parameters (generator matrices, etc.)
    """

    # Mapping from distribution type strings to sampler classes
    DISTRIBUTION_TYPES = {
        'EXP': ExponentialSampler,
        'EXPONENTIAL': ExponentialSampler,
        'ERLANG': ErlangSampler,
        'GAMMA': GammaSampler,
        'UNIFORM': UniformSampler,
        'DET': DeterministicSampler,
        'DETERMINISTIC': DeterministicSampler,
        'PARETO': ParetoSampler,
        'WEIBULL': WeibullSampler,
        'LOGNORMAL': LognormalSampler,
        'NORMAL': NormalSampler,
        'BETA': BetaSampler,
        'PH': PHSampler,
        'APH': APHSampler,
        'COXIAN': CoxianSampler,
        'COX2': Coxian2Sampler,
        'HYPEREXP': HyperExponentialSampler,
        'HYPEREXP2': HyperExp2Sampler,
        'MAP': MAPSampler,
        'MAP2': MAP2Sampler,
        'BMAP': BMAPSampler,
        'MMPP': MMPPSampler,
        'MMPP2': MMPPSampler,
        'IMMEDIATE': ImmediateSampler,
        'DISABLED': DisabledSampler,
        # Discrete distributions
        'POISSON': PoissonSampler,
        'BINOMIAL': BinomialSampler,
        'GEOMETRIC': GeometricSampler,
        'NEGBINOMIAL': NegativeBinomialSampler,
        'NEGATIVEBINOMIAL': NegativeBinomialSampler,
        'BERNOULLI': BernoulliSampler,
        'DISCRETEUNIFORM': DiscreteUniformSampler,
        'DUNIFORM': DiscreteUniformSampler,
        'ZIPF': ZipfSampler,
        'CATEGORICAL': CategoricalSampler,
        'EMPIRICAL': EmpiricalDiscreteSampler,
        'EMPIRICALDISCRETE': EmpiricalDiscreteSampler,
        # Trace replay
        'TRACE': FileTraceSampler,
        'FILETRACE': FileTraceSampler,
        'REPLAY': FileTraceSampler,
    }

    def __init__(self, rng: Optional[np.random.Generator] = None):
        """
        Initialize the distribution factory.

        Args:
            rng: Random number generator (optional)
        """
        self.rng = rng if rng is not None else np.random.default_rng()

    def create_from_rate(self, rate: float) -> DistributionSampler:
        """
        Create exponential sampler from rate parameter.

        For backward compatibility, rate-based specifications are
        interpreted as exponential distributions.

        Args:
            rate: Rate parameter (must be positive)

        Returns:
            ExponentialSampler configured with the given rate

        Raises:
            ValueError: If rate is not positive
        """
        if rate <= 0:
            raise ValueError(f"Rate must be positive, got {rate}")
        return ExponentialSampler(rate=rate, rng=self.rng)

    def create_from_spec(
        self,
        dist_type: str,
        parameters: Dict[str, Any]
    ) -> DistributionSampler:
        """
        Create distribution sampler from explicit type and parameters.

        Args:
            dist_type: Distribution type string (EXP, ERLANG, PH, MAP, etc.)
            parameters: Dictionary of distribution parameters

        Returns:
            Configured distribution sampler

        Raises:
            ValueError: If distribution type is unknown or parameters invalid
            KeyError: If required parameters are missing
        """
        dist_type_upper = dist_type.upper().strip()

        if dist_type_upper not in self.DISTRIBUTION_TYPES:
            raise ValueError(
                f"Unknown distribution type: {dist_type}. "
                f"Supported types: {', '.join(sorted(self.DISTRIBUTION_TYPES.keys()))}"
            )

        sampler_class = self.DISTRIBUTION_TYPES[dist_type_upper]

        # Handle special cases
        if dist_type_upper == 'IMMEDIATE':
            return ImmediateSampler(rng=self.rng)

        if dist_type_upper == 'DISABLED':
            return DisabledSampler(rng=self.rng)

        # Handle exponential (single rate parameter)
        if dist_type_upper in ['EXP', 'EXPONENTIAL']:
            if 'rate' in parameters:
                return ExponentialSampler(rate=parameters['rate'], rng=self.rng)
            elif 'mean' in parameters:
                return ExponentialSampler(rate=1.0 / parameters['mean'], rng=self.rng)
            else:
                raise KeyError(f"EXP requires 'rate' or 'mean' parameter")

        # Handle Erlang (shape, rate)
        if dist_type_upper == 'ERLANG':
            if 'shape' not in parameters or 'rate' not in parameters:
                raise KeyError("ERLANG requires 'shape' and 'rate' parameters")
            return ErlangSampler(
                shape=int(parameters['shape']),
                rate=parameters['rate'],
                rng=self.rng
            )

        # Handle Gamma (shape, rate)
        if dist_type_upper == 'GAMMA':
            if 'shape' not in parameters or 'rate' not in parameters:
                raise KeyError("GAMMA requires 'shape' and 'rate' parameters")
            return GammaSampler(
                shape=parameters['shape'],
                rate=parameters['rate'],
                rng=self.rng
            )

        # Handle Uniform (min, max)
        if dist_type_upper == 'UNIFORM':
            if 'min' not in parameters or 'max' not in parameters:
                raise KeyError("UNIFORM requires 'min' and 'max' parameters")
            return UniformSampler(
                min_val=parameters['min'],
                max_val=parameters['max'],
                rng=self.rng
            )

        # Handle Deterministic (value)
        if dist_type_upper in ['DET', 'DETERMINISTIC']:
            if 'value' not in parameters:
                raise KeyError("DETERMINISTIC requires 'value' parameter")
            return DeterministicSampler(value=parameters['value'], rng=self.rng)

        # Handle Pareto (scale, shape)
        if dist_type_upper == 'PARETO':
            if 'scale' not in parameters or 'shape' not in parameters:
                raise KeyError("PARETO requires 'scale' and 'shape' parameters")
            return ParetoSampler(
                scale=parameters['scale'],
                shape=parameters['shape'],
                rng=self.rng
            )

        # Handle Weibull (shape, scale)
        if dist_type_upper == 'WEIBULL':
            if 'shape' not in parameters or 'scale' not in parameters:
                raise KeyError("WEIBULL requires 'shape' and 'scale' parameters")
            return WeibullSampler(
                shape=parameters['shape'],
                scale=parameters['scale'],
                rng=self.rng
            )

        # Handle Lognormal (mu, sigma)
        if dist_type_upper == 'LOGNORMAL':
            if 'mu' not in parameters or 'sigma' not in parameters:
                raise KeyError("LOGNORMAL requires 'mu' and 'sigma' parameters")
            return LognormalSampler(
                mu=parameters['mu'],
                sigma=parameters['sigma'],
                rng=self.rng
            )

        # Handle Normal (mu, sigma, min, max for truncation)
        if dist_type_upper == 'NORMAL':
            if 'mu' not in parameters or 'sigma' not in parameters:
                raise KeyError("NORMAL requires 'mu' and 'sigma' parameters")
            return NormalSampler(
                mu=parameters['mu'],
                sigma=parameters['sigma'],
                min_val=parameters.get('min'),
                max_val=parameters.get('max'),
                rng=self.rng
            )

        # Handle Beta (alpha, beta)
        if dist_type_upper == 'BETA':
            if 'alpha' not in parameters or 'beta' not in parameters:
                raise KeyError("BETA requires 'alpha' and 'beta' parameters")
            return BetaSampler(
                alpha=parameters['alpha'],
                beta=parameters['beta'],
                rng=self.rng
            )

        # Handle Phase-Type (moments or matrix)
        if dist_type_upper == 'PH':
            if 'moments' in parameters:
                return PHSampler.from_moments(
                    moments=parameters['moments'],
                    rng=self.rng
                )
            elif 'matrix' in parameters:
                return PHSampler(
                    T=parameters['matrix'],
                    rng=self.rng
                )
            else:
                raise KeyError("PH requires 'moments' or 'matrix' parameter")

        # Handle Acyclic Phase-Type (moments or matrix)
        if dist_type_upper == 'APH':
            if 'moments' in parameters:
                return APHSampler.from_moments(
                    moments=parameters['moments'],
                    rng=self.rng
                )
            elif 'matrix' in parameters:
                return APHSampler(
                    T=parameters['matrix'],
                    rng=self.rng
                )
            else:
                raise KeyError("APH requires 'moments' or 'matrix' parameter")

        # Handle Coxian (moments)
        if dist_type_upper == 'COXIAN':
            if 'moments' in parameters:
                return CoxianSampler.from_moments(
                    moments=parameters['moments'],
                    rng=self.rng
                )
            else:
                raise KeyError("COXIAN requires 'moments' parameter")

        # Handle 2-stage Coxian (mean, scv)
        if dist_type_upper == 'COX2':
            if 'mean' not in parameters or 'scv' not in parameters:
                raise KeyError("COX2 requires 'mean' and 'scv' parameters")
            return Coxian2Sampler(
                mean=parameters['mean'],
                scv=parameters['scv'],
                rng=self.rng
            )

        # Handle HyperExponential
        if dist_type_upper == 'HYPEREXP':
            if 'moments' in parameters:
                return HyperExponentialSampler.from_moments(
                    moments=parameters['moments'],
                    rng=self.rng
                )
            else:
                raise KeyError("HYPEREXP requires 'moments' parameter")

        # Handle 2-stage HyperExponential (mean, scv)
        if dist_type_upper == 'HYPEREXP2':
            if 'mean' not in parameters or 'scv' not in parameters:
                raise KeyError("HYPEREXP2 requires 'mean' and 'scv' parameters")
            return HyperExp2Sampler(
                mean=parameters['mean'],
                scv=parameters['scv'],
                rng=self.rng
            )

        # Handle MAP/BMAP (generator matrices or parameters)
        if dist_type_upper == 'MAP':
            if 'D0' in parameters and 'D1' in parameters:
                return MAPSampler(
                    D0=parameters['D0'],
                    D1=parameters['D1'],
                    rng=self.rng
                )
            else:
                raise KeyError("MAP requires 'D0' and 'D1' matrices")

        if dist_type_upper in ['MAP2', 'MMPP', 'MMPP2']:
            if 'D0' in parameters and 'D1' in parameters:
                return MAP2Sampler(
                    D0=parameters['D0'],
                    D1=parameters['D1'],
                    rng=self.rng
                )
            else:
                raise KeyError(f"{dist_type_upper} requires 'D0' and 'D1' matrices")

        if dist_type_upper == 'BMAP':
            if 'D0' in parameters and 'D' in parameters:
                return BMAPSampler(
                    D0=parameters['D0'],
                    D=parameters['D'],
                    rng=self.rng
                )
            else:
                raise KeyError("BMAP requires 'D0' and 'D' matrices")

        # Handle Poisson (lambda)
        if dist_type_upper == 'POISSON':
            if 'lambda' in parameters:
                return PoissonSampler(lambda_=parameters['lambda'], rng=self.rng)
            elif 'rate' in parameters:
                return PoissonSampler(lambda_=parameters['rate'], rng=self.rng)
            elif 'mean' in parameters:
                return PoissonSampler(lambda_=parameters['mean'], rng=self.rng)
            else:
                raise KeyError("POISSON requires 'lambda', 'rate', or 'mean' parameter")

        # Handle Binomial (n, p)
        if dist_type_upper == 'BINOMIAL':
            if 'n' not in parameters or 'p' not in parameters:
                raise KeyError("BINOMIAL requires 'n' and 'p' parameters")
            return BinomialSampler(
                n=int(parameters['n']),
                p=parameters['p'],
                rng=self.rng
            )

        # Handle Geometric (p)
        if dist_type_upper == 'GEOMETRIC':
            if 'p' not in parameters:
                raise KeyError("GEOMETRIC requires 'p' parameter")
            shift = parameters.get('shift', True)
            return GeometricSampler(p=parameters['p'], shift=shift, rng=self.rng)

        # Handle Negative Binomial (r, p)
        if dist_type_upper in ['NEGBINOMIAL', 'NEGATIVEBINOMIAL']:
            if 'r' not in parameters or 'p' not in parameters:
                raise KeyError("NEGATIVEBINOMIAL requires 'r' and 'p' parameters")
            return NegativeBinomialSampler(
                r=parameters['r'],
                p=parameters['p'],
                rng=self.rng
            )

        # Handle Bernoulli (p)
        if dist_type_upper == 'BERNOULLI':
            if 'p' not in parameters:
                raise KeyError("BERNOULLI requires 'p' parameter")
            return BernoulliSampler(p=parameters['p'], rng=self.rng)

        # Handle Discrete Uniform (low, high)
        if dist_type_upper in ['DISCRETEUNIFORM', 'DUNIFORM']:
            if 'low' not in parameters or 'high' not in parameters:
                raise KeyError("DISCRETEUNIFORM requires 'low' and 'high' parameters")
            return DiscreteUniformSampler(
                low=int(parameters['low']),
                high=int(parameters['high']),
                rng=self.rng
            )

        # Handle Zipf (s, n)
        if dist_type_upper == 'ZIPF':
            if 's' not in parameters:
                raise KeyError("ZIPF requires 's' (exponent) parameter")
            return ZipfSampler(
                s=parameters['s'],
                n=parameters.get('n'),
                rng=self.rng
            )

        # Handle Categorical (probabilities, values)
        if dist_type_upper == 'CATEGORICAL':
            if 'probabilities' not in parameters:
                raise KeyError("CATEGORICAL requires 'probabilities' parameter")
            return CategoricalSampler(
                probabilities=parameters['probabilities'],
                values=parameters.get('values'),
                rng=self.rng
            )

        # Handle Empirical Discrete (values, counts)
        if dist_type_upper in ['EMPIRICAL', 'EMPIRICALDISCRETE']:
            if 'values' not in parameters:
                raise KeyError("EMPIRICALDISCRETE requires 'values' parameter")
            return EmpiricalDiscreteSampler(
                values=parameters['values'],
                counts=parameters.get('counts'),
                rng=self.rng
            )

        # Handle Trace File (file_path, column, loop)
        if dist_type_upper in ['TRACE', 'FILETRACE', 'REPLAY']:
            if 'file_path' not in parameters and 'path' not in parameters:
                raise KeyError("TRACE requires 'file_path' or 'path' parameter")
            file_path = parameters.get('file_path', parameters.get('path'))
            return FileTraceSampler(
                file_path=file_path,
                column=parameters.get('column', 0),
                skip_header=parameters.get('skip_header', False),
                loop=parameters.get('loop', True),
                delimiter=parameters.get('delimiter'),
                rng=self.rng
            )

        # Fallback (should not reach here)
        raise NotImplementedError(
            f"Distribution factory does not support {dist_type_upper} yet"
        )

    def create_arrival_samplers(
        self,
        sn: Any
    ) -> Dict[Tuple[int, int], DistributionSampler]:
        """
        Create arrival distribution samplers from NetworkStruct.

        Extracts arrival distribution specifications for each source node
        and job class, returning a dictionary indexed by (source_idx, class_id).

        Args:
            sn: NetworkStruct with arrival distribution specifications

        Returns:
            Dictionary: (source_idx, class_id) -> DistributionSampler

        Notes:
            - Requires: sn.arrivals, sn.arrivalProcesses (if defined)
            - Falls back to rate-based exponential if explicit distributions unavailable
        """
        samplers = {}

        # First try to get explicit arrival process distributions
        if hasattr(sn, 'arrivalProcesses') and sn.arrivalProcesses is not None:
            try:
                arrival_processes = sn.arrivalProcesses
                num_sources = len([n for n in range(len(sn.nodetypes))
                                  if sn.nodetypes[n] == 0])  # SOURCE = 0
                num_classes = sn.nclasses if hasattr(sn, 'nclasses') else 1

                for src_idx in range(num_sources):
                    for class_id in range(num_classes):
                        # Try to get distribution spec for this (source, class)
                        try:
                            dist_spec = arrival_processes.get((src_idx, class_id))
                            if dist_spec is not None and hasattr(dist_spec, 'type'):
                                dist_type = str(dist_spec.type).upper()
                                params = dist_spec.params if hasattr(dist_spec, 'params') else {}
                                samplers[(src_idx, class_id)] = self.create_from_spec(dist_type, params)
                        except Exception:
                            pass  # Skip on error, will use rate-based fallback
            except Exception:
                pass  # If arrivalProcesses can't be parsed, use rates

        # Otherwise, create exponential samplers from arrival rates (fallback)
        if hasattr(sn, 'arrivals') and sn.arrivals is not None:
            arrivals = sn.arrivals
            num_sources = len([n for n in range(len(sn.nodetypes))
                              if sn.nodetypes[n] == 0])  # SOURCE = 0
            num_classes = sn.nclasses

            for src_idx in range(num_sources):
                for class_id in range(num_classes):
                    # Only create if not already in samplers
                    if (src_idx, class_id) not in samplers and arrivals[src_idx, class_id] > 0:
                        samplers[(src_idx, class_id)] = self.create_from_rate(
                            arrivals[src_idx, class_id]
                        )

        return samplers

    def create_service_samplers(
        self,
        sn: Any
    ) -> Dict[Tuple[int, int], DistributionSampler]:
        """
        Create service distribution samplers from NetworkStruct.

        Extracts service distribution specifications for each service node
        and job class, returning a dictionary indexed by (station_idx, class_id).

        Args:
            sn: NetworkStruct with service distribution specifications

        Returns:
            Dictionary: (station_idx, class_id) -> DistributionSampler

        Notes:
            - Requires: sn.rates, sn.proc, sn.procid
            - Falls back to rate-based exponential if explicit distributions unavailable
        """
        samplers = {}
        num_stations = sn.nstations
        num_classes = sn.nclasses

        # Check for explicit process distributions (proc/procid fields)
        has_proc = hasattr(sn, 'proc') and sn.proc is not None
        has_procid = hasattr(sn, 'procid') and sn.procid is not None

        for stn_idx in range(num_stations):
            for class_id in range(num_classes):
                sampler = None

                # Try to create sampler from proc/procid
                if has_proc and has_procid:
                    sampler = self._create_sampler_from_proc(sn, stn_idx, class_id)

                # Fallback to rate-based exponential
                if sampler is None:
                    if hasattr(sn, 'rates') and sn.rates is not None:
                        if stn_idx < sn.rates.shape[0] and class_id < sn.rates.shape[1]:
                            rate = sn.rates[stn_idx, class_id]
                            if rate > 0 and not np.isnan(rate):
                                sampler = self.create_from_rate(rate)

                if sampler is not None:
                    samplers[(stn_idx, class_id)] = sampler

        return samplers

    def _create_sampler_from_proc(
        self,
        sn: Any,
        stn_idx: int,
        class_id: int
    ) -> Optional[DistributionSampler]:
        """
        Create a distribution sampler from proc/procid fields.

        Args:
            sn: NetworkStruct with proc and procid
            stn_idx: Station index
            class_id: Class index

        Returns:
            DistributionSampler or None if not available
        """
        try:
            # Get process type
            procid = None
            if hasattr(sn, 'procid') and sn.procid is not None:
                if stn_idx < sn.procid.shape[0] and class_id < sn.procid.shape[1]:
                    procid = sn.procid[stn_idx, class_id]

            if procid is None:
                return None

            # Get process parameters
            proc = None
            if hasattr(sn, 'proc') and sn.proc is not None:
                try:
                    if stn_idx < len(sn.proc) and sn.proc[stn_idx] is not None:
                        if class_id < len(sn.proc[stn_idx]):
                            proc = sn.proc[stn_idx][class_id]
                except (IndexError, TypeError):
                    pass

            # Import ProcessType for comparison
            from ....constants import ProcessType

            # Handle MAP/MMPP2
            if procid in (ProcessType.MAP, ProcessType.MMPP2):
                if proc is not None and isinstance(proc, (list, tuple)) and len(proc) >= 2:
                    D0 = np.asarray(proc[0], dtype=np.float64)
                    D1 = np.asarray(proc[1], dtype=np.float64)
                    return MAP2Sampler(D0=D0, D1=D1, rng=self.rng)

            # Handle PH/APH
            if procid in (ProcessType.PH, ProcessType.APH):
                if proc is not None and isinstance(proc, (list, tuple)) and len(proc) >= 2:
                    alpha = np.asarray(proc[0], dtype=np.float64)
                    T = np.asarray(proc[1], dtype=np.float64)
                    return PHSampler(T=T, alpha=alpha, rng=self.rng)

            # Handle COXIAN
            if procid == ProcessType.COXIAN:
                if proc is not None and isinstance(proc, (list, tuple)) and len(proc) >= 2:
                    alpha = np.asarray(proc[0], dtype=np.float64)
                    T = np.asarray(proc[1], dtype=np.float64)
                    return CoxianSampler(T=T, alpha=alpha, rng=self.rng)

            # Handle HYPEREXP
            if procid == ProcessType.HYPEREXP:
                if proc is not None and isinstance(proc, dict):
                    p = proc.get('p', 0.5)
                    mu1 = proc.get('mu1', 1.0)
                    mu2 = proc.get('mu2', 1.0)
                    return HyperExp2Sampler(p=p, mu1=mu1, mu2=mu2, rng=self.rng)

            # Handle ERLANG
            if procid == ProcessType.ERLANG:
                if proc is not None and isinstance(proc, dict):
                    k = int(proc.get('k', 1))
                    mu = proc.get('mu', 1.0)
                    return ErlangSampler(shape=k, rate=mu, rng=self.rng)

            # Handle EXP (default)
            if procid == ProcessType.EXP:
                if proc is not None and isinstance(proc, dict):
                    rate = proc.get('rate', 1.0)
                    return ExponentialSampler(rate=rate, rng=self.rng)

            # For other types, return None to fall back to rate-based
            return None

        except Exception:
            # On any error, return None to fall back to rate-based
            return None

    def create_iats_samplers(
        self,
        sn: Any
    ) -> Dict[Tuple[int, int], DistributionSampler]:
        """
        Create inter-arrival time distribution samplers from NetworkStruct.

        Similar to create_arrival_samplers but works with IAT (inter-arrival times)
        instead of rates.

        Args:
            sn: NetworkStruct with IAT specifications

        Returns:
            Dictionary: (source_idx, class_id) -> DistributionSampler
        """
        samplers = {}

        # If IAT matrix available, use it
        if hasattr(sn, 'iats') and sn.iats is not None:
            iats = sn.iats
            num_sources = len([n for n in range(len(sn.nodetypes))
                              if sn.nodetypes[n] == 0])  # SOURCE = 0
            num_classes = sn.nclasses

            for src_idx in range(num_sources):
                for class_id in range(num_classes):
                    if iats[src_idx, class_id] > 0:
                        # IAT to rate: rate = 1 / IAT
                        rate = 1.0 / iats[src_idx, class_id]
                        samplers[(src_idx, class_id)] = self.create_from_rate(rate)

        return samplers

    def create_service_times_samplers(
        self,
        sn: Any
    ) -> Dict[Tuple[int, int], DistributionSampler]:
        """
        Create service time distribution samplers from NetworkStruct.

        Similar to create_service_samplers but works with service times
        instead of rates.

        Args:
            sn: NetworkStruct with service time specifications

        Returns:
            Dictionary: (station_idx, class_id) -> DistributionSampler
        """
        samplers = {}

        # If service times matrix available, use it
        if hasattr(sn, 'serviceTimes') and sn.serviceTimes is not None:
            service_times = sn.serviceTimes
            num_stations = sn.nstations
            num_classes = sn.nclasses

            for stn_idx in range(num_stations):
                for class_id in range(num_classes):
                    if service_times[stn_idx, class_id] > 0:
                        # Service time to rate: rate = 1 / mean_service_time
                        rate = 1.0 / service_times[stn_idx, class_id]
                        samplers[(stn_idx, class_id)] = self.create_from_rate(rate)

        return samplers
