
import sys

import jpype
import jpype.imports
import numpy as np

from line_solver import jlineMatrixFromArray, jlineMatrixToArray

class NamedParam:
    """
    Named parameter for probability distributions.
    
    This class represents a parameter with both a name and value,
    used to configure probability distributions.
    
    Attributes:
        name (str): Name of the parameter.
        value: Value of the parameter.
    """
    
    def __init__(self, *args):
        if len(args) == 1:
            self.obj = args[0]
        else:
            self.name = args[0]
            self.value = args[1]

    def getName(self):
        """
        Get the name of this named parameter.
        
        Returns:
            str: The parameter name.
        """
        return self.name

    def getValue(self):
        """
        Get the value of this named parameter.
        
        Returns:
            The parameter value (type depends on the specific parameter).
        """
        return self.value

    get_name = getName
    get_value = getValue

class Distribution:
    """
    Base class for all probability distributions.
    
    This class provides the common interface for all probability distributions
    in LINE, including service time distributions, inter-arrival time distributions,
    and other stochastic processes.
    
    The class supports:
    - Continuous and discrete distributions
    - Markovian and non-Markovian distributions
    - Statistical moments (mean, variance, skewness)
    - Cumulative distribution function (CDF)
    - Laplace-Stieltjes transform (LST)
    - Random sampling
    """

    def __init__(self):
        """Initialize a new distribution."""
        pass

    def evalCDF(self, x):
        """
        Evaluate the cumulative distribution function at point x.
        
        Args:
            x (float): Point at which to evaluate the CDF.
            
        Returns:
            float: Value of F(x) = P(X <= x).
        """
        return self.obj.evalCDF(x)

    def evalLST(self, x):
        """
        Evaluate the Laplace-Stieltjes transform at point x.
        
        Args:
            x (float): Point at which to evaluate the LST.
            
        Returns:
            float: Value of the Laplace-Stieltjes transform.
        """
        return self.obj.evalLST(x)

    def getName(self):
        """
        Get the name of this distribution.
        
        Returns:
            str: Distribution name (e.g., 'Exp', 'Det', 'Erlang').
        """
        return self.obj.getName()

    def getParam(self, id):
        """
        Get a named parameter of this distribution by its ID.
        
        Args:
            id (int): The parameter ID (0-based index).
            
        Returns:
            NamedParam: A named parameter object containing the parameter name and value.
        """
        nparam = NamedParam(self.obj.getParam(id))
        return nparam

    def getMean(self):
        """
        Get the mean (expected value) of the distribution.
        
        Returns:
            float: Mean value E[X].
        """
        return self.obj.getMean()

    def getRate(self):
        """
        Get the rate parameter (1/mean) of the distribution.
        
        Returns:
            float: Rate parameter λ = 1/E[X].
        """
        return self.obj.getRate()

    def getSCV(self):
        """
        Get the squared coefficient of variation (SCV).
        
        The SCV is the ratio of variance to squared mean: Var[X]/E[X]².
        
        Returns:
            float: Squared coefficient of variation.
        """
        return self.obj.getSCV()

    def getVar(self):
        """
        Get the variance of the distribution.
        
        Returns:
            float: Variance Var[X] = E[X²] - E[X]².
        """
        return self.obj.getVar()

    def getSkew(self):
        """
        Get the skewness of the distribution.
        
        Skewness measures the asymmetry of the distribution.
        
        Returns:
            float: Skewness coefficient (0 for symmetric distributions).
        """
        return self.obj.getSkew()

    def getSupport(self):
        """
        Get the support of the distribution.
        
        The support is the set of values for which the distribution 
        has non-zero probability density/mass.
        
        Returns:
            Matrix: A 2-element array [min, max] defining the support range.
        """
        return self.obj.getSupport()

    def isContinuous(self):
        """
        Check if this distribution is continuous.
        
        Returns:
            bool: True if the distribution is continuous, False otherwise.
        """
        return self.obj.isContinuous()

    def isDisabled(self):
        """
        Check if this distribution is disabled.
        
        A disabled distribution represents a null or inactive service.
        
        Returns:
            bool: True if the distribution is disabled, False otherwise.
        """
        return self.obj.isDisabled()

    def isDiscrete(self):
        """
        Check if this distribution is discrete.
        
        Returns:
            bool: True if the distribution is discrete, False otherwise.
        """
        return self.obj.isDiscrete()

    def isImmediate(self):
        """
        Check if this distribution represents immediate service.
        
        An immediate distribution has zero service time.
        
        Returns:
            bool: True if service is immediate, False otherwise.
        """
        return self.obj.isImmediate()

    def sample(self, *args):
        """
        Generate random samples from this distribution.
        
        Args:
            n (int): Number of samples to generate.
            seed (int, optional): Random seed for reproducible sampling.
            
        Returns:
            list: Array of n random samples from the distribution.
        """
        if len(args) == 1:
            n = args[0]
            return jlineMatrixToArray(self.obj.sample(n))
        else:
            n = args[0]
            seed = args[1]

    get_name = getName
    get_param = getParam
    get_mean = getMean
    get_rate = getRate
    get_scv = getSCV
    get_var = getVar
    get_skew = getSkew
    get_support = getSupport

class ContinuousDistribution(Distribution):
    """
    Base class for continuous probability distributions.
    
    Continuous distributions have support over real number intervals
    and are commonly used for service times and inter-arrival times
    in queueing models.
    """
    
    def __init__(self):
        """Initialize a continuous distribution."""
        super().__init__()

class DiscreteDistribution(Distribution):
    """
    Base class for discrete probability distributions.
    
    Discrete distributions have support over discrete values
    (integers or specific points) and are used for batch sizes,
    discrete service requirements, and other countable quantities.
    """
    
    def __init__(self):
        """Initialize a discrete distribution."""
        super().__init__()

class Markovian(Distribution):
    def __init__(self):
        super().__init__()

    def getD0(self):
        return self.obj.getD0()

    def getD1(self):
        return self.obj.getD1()

    def getMu(self):
        return self.obj.getMu()

    def getNumberOfPhases(self):
        return self.obj.getNumberOfPhases()

    def getPH(self):
        return self.obj.getPH()

    def getPhi(self):
        return self.obj.getPhi()

    def getInitProb(self):
        return self.obj.getInitProb()

    get_d0 = getD0
    get_d1 = getD1
    get_mu = getMu
    get_number_of_phases = getNumberOfPhases
    get_ph = getPH
    get_phi = getPhi
    get_init_prob = getInitProb

class APH(Markovian):
    """
    Acyclic Phase-type (APH) distribution.
    
    Represents an acyclic phase-type distribution defined by initial
    probability vector and sub-generator matrix. APH distributions
    are a subclass of PH distributions with acyclic structure.
    
    Args:
        alpha: Initial probability vector or existing APH object
        subgen: Sub-generator matrix (if alpha is vector)
    """
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            self.obj = args[0]
        else:
            alpha = args[0]
            subgen = args[1]
            self.obj = jpype.JPackage('jline').lang.processes.APH(jlineMatrixFromArray(alpha),
                                                                      jlineMatrixFromArray(subgen))

    @staticmethod
    def fitMeanAndSCV(mean, scv):
        return APH(jpype.JPackage('jline').lang.processes.APH.fitMeanAndSCV(mean, scv))

    fit_mean_and_scv = fitMeanAndSCV

    @staticmethod
    def fitCentral(mean, var, skew):
        """
        Create an APH distribution matching the first three central moments.

        Args:
            mean (float): Mean (first moment).
            var (float): Variance (second central moment).
            skew (float): Skewness (third standardized moment).

        Returns:
            APH: An APH distribution matching the given moments.
        """
        return APH(jpype.JPackage('jline').lang.processes.APH.fitCentral(mean, var, skew))

    fit_central = fitCentral

class Bernoulli(DiscreteDistribution):
    """
    Bernoulli distribution for binary outcomes.
    
    Represents a discrete distribution that takes value 1 with probability p
    and value 0 with probability (1-p).
    
    Args:
        prob (float): Success probability p, must be in [0,1].
    """
    def __init__(self, *args):
        super().__init__()
        prob = args[0]
        self.obj = jpype.JPackage('jline').lang.processes.Bernoulli(prob)

class Binomial(DiscreteDistribution):
    """
    Binomial distribution for number of successes in n trials.
    
    Represents the distribution of the number of successes in n independent
    Bernoulli trials, each with success probability p.
    
    Args:
        prob (float): Success probability for each trial.
        n (int): Number of trials.
    """
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            self.obj = args[0]
        else:
            prob = args[0]
            n = args[1]
            self.obj = jpype.JPackage('jline').lang.processes.Binomial(prob, n)

class Coxian(Markovian):
    """
    Coxian distribution for phase-type service times.
    
    A Coxian distribution is a special case of phase-type distribution
    with a series structure, commonly used for modeling service times
    with controlled variance.
    
    Args:
        mu (array_like): Service rates for each phase.
        phi (array_like): Exit probabilities from each phase.
    """
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            self.obj = args[0]
        else:
            mu = args[0]
            phi = args[1]

            if hasattr(phi, 'toArray'):
                phi_array = phi.toArray()
            else:
                phi_array = phi

            if hasattr(mu, 'obj'):
                mu_matrix = mu.obj
            else:
                mu_matrix = jlineMatrixFromArray(mu)

            if hasattr(phi, 'obj'):
                phi_matrix = phi.obj
            else:
                phi_matrix = jlineMatrixFromArray(phi)

            if len(phi_array.shape) == 2:
                phi_flat = phi_array.flatten()
            else:
                phi_flat = phi_array

            if phi_flat[-1] != 1.0:
                print("Invalid Coxian exit probabilities. The last element must be 1.0.", file=sys.stderr)
            elif max(phi_flat) > 1.0 or min(phi_flat) < 0.0:
                print("Invalid Coxian exit probabilities. Some values are not in [0,1].", file=sys.stderr)
            else:
                self.obj = jpype.JPackage('jline').lang.processes.Coxian(mu_matrix, phi_matrix)

    @staticmethod
    def fitMeanAndSCV(mean, scv):
        """
        Create a Coxian distribution with specified mean and squared coefficient of variation.

        Args:
            mean (float): Mean service time > 0.
            scv (float): Squared coefficient of variation.

        Returns:
            Coxian: A Coxian distribution matching the given mean and SCV.
        """
        return Coxian(jpype.JPackage('jline').lang.processes.Coxian.fitMeanAndSCV(mean, scv))

    fit_mean_and_scv = fitMeanAndSCV

    @staticmethod
    def fitCentral(mean, var, skew):
        """
        Create a Coxian distribution matching the first three central moments.

        Args:
            mean (float): Mean (first moment).
            var (float): Variance (second central moment).
            skew (float): Skewness (third standardized moment).

        Returns:
            Coxian: A Coxian distribution matching the given moments.
        """
        return Coxian(jpype.JPackage('jline').lang.processes.Coxian.fitCentral(mean, var, skew))

    fit_central = fitCentral

class Cox2(Markovian):
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            self.obj = args[0]
        else:
            mu1 = args[0]
            mu2 = args[1]
            phi1 = args[2]
            self.obj = jpype.JPackage('jline').lang.processes.Cox2(mu1, mu2, phi1)

    @staticmethod
    def fitMeanAndSCV(mean, scv):
        """
        Create a Cox2 distribution with specified mean and squared coefficient of variation.

        Args:
            mean (float): Mean service time > 0.
            scv (float): Squared coefficient of variation.

        Returns:
            Cox2: A Cox2 distribution matching the given mean and SCV.
        """
        return Cox2(jpype.JPackage('jline').lang.processes.Cox2.fitMeanAndSCV(mean, scv))

    fit_mean_and_scv = fitMeanAndSCV

    @staticmethod
    def fitCentral(mean, var, skew):
        """
        Create a Cox2 distribution matching the first three central moments.

        Args:
            mean (float): Mean (first moment).
            var (float): Variance (second central moment).
            skew (float): Skewness (third standardized moment).

        Returns:
            Cox2: A Cox2 distribution matching the given moments.
        """
        return Cox2(jpype.JPackage('jline').lang.processes.Cox2.fitCentral(mean, var, skew))

    fit_central = fitCentral


class Det(Distribution):
    """
    Deterministic (constant) distribution.
    
    The deterministic distribution represents a constant value with no
    randomness. All samples return the same fixed value. It has zero
    variance and SCV = 0.
    
    This distribution is often used for:
    - Fixed service times
    - Deterministic delays
    - Constant inter-arrival times in fluid models
    
    Args:
        value (float): The constant value.
        
    Examples:
        >>> det_dist = Det(5.0)  # Constant value = 5.0
        >>> mean = det_dist.getMean()  # Returns 5.0
        >>> var = det_dist.getVar()    # Returns 0.0
        >>> scv = det_dist.getSCV()    # Returns 0.0
    """

    def __init__(self, value):
        super().__init__()
        if isinstance(value, (int, float, np.integer, np.floating)):
            self.obj = jpype.JPackage('jline').lang.processes.Det(float(value))
        else:
            self.obj = value

    @staticmethod
    def fitMean(mean):
        """
        Create a deterministic distribution with the specified mean.

        Since Det has zero variance, the mean equals the constant value.

        Args:
            mean (float): The constant value for the distribution.

        Returns:
            Det: A deterministic distribution with value = mean.
        """
        return Det(jpype.JPackage('jline').lang.processes.Det.fitMean(mean))

    fit_mean = fitMean


class Disabled(Distribution):
    def __init__(self, *args):
        super().__init__()
        if len(args) == 0:
            self.obj = jpype.JPackage('jline').lang.processes.Disabled()
        else:
            self.obj = args[0]

    @staticmethod
    def getInstance():
        return Disabled(jpype.JPackage('jline').lang.processes.Disabled.getInstance())

    get_instance = getInstance

class DiscreteSampler(DiscreteDistribution):
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            if isinstance(args[0], DiscreteDistribution):
                self.obj = args[0]
            else:
                p = args[0]
                if isinstance(p, list):
                    p = jlineMatrixFromArray(p)
                elif hasattr(p, 'obj'):
                    p = p.obj
                self.obj = jpype.JPackage('jline').lang.processes.DiscreteSampler(p)
        else:
            p = args[0]
            x = args[1]
            if isinstance(p, list):
                p = jlineMatrixFromArray(p)
            elif hasattr(p, 'obj'):
                p = p.obj
            if isinstance(x, list):
                x = jlineMatrixFromArray(x)
            elif hasattr(x, 'obj'):
                x = x.obj
            self.obj = jpype.JPackage('jline').lang.processes.DiscreteSampler(p, x)

class DiscreteUniform(DiscreteDistribution):
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            self.obj = args[0]
        else:
            minVal = args[0]
            maxVal = args[1]
            self.obj = jpype.JPackage('jline').lang.processes.DiscreteUniform(minVal, maxVal)

class Exp(Markovian):
    """
    Exponential distribution for memoryless service and arrival processes.
    
    The exponential distribution is the most common distribution used in queueing
    theory due to its memoryless property. It's characterized by a single rate parameter.
    
    Args:
        rate (float): Rate parameter (λ). Mean service time = 1/rate.
        
    Examples:
        >>> Exp(2.0)  # Exponential with rate 2.0, mean = 0.5
        >>> Exp.fitMean(0.5)  # Exponential fitted to mean 0.5
    """

    def __init__(self, *args):
        super().__init__()

        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, (int, float, np.integer, np.floating)):
                self.obj = jpype.JPackage('jline').lang.processes.Exp.fitRate(arg)
            else:
                self.obj = arg
        else:
            raise ValueError("Exp constructor accepts a single rate (float) or a pre-constructed object.")

    def fitRate(rate):
        """
        Create an exponential distribution with the specified rate parameter.
        
        Args:
            rate (float): Rate parameter λ > 0. Higher values mean faster service.
            
        Returns:
            Exp: An exponential distribution with mean = 1/rate.
        """
        return Exp(jpype.JPackage('jline').lang.processes.Exp.fitRate(rate))

    fit_rate = fitRate

    def fitMean(mean):
        """
        Create an exponential distribution with the specified mean.
        
        Args:
            mean (float): Mean service time > 0. The rate will be 1/mean.
            
        Returns:
            Exp: An exponential distribution with the specified mean.
        """
        return Exp(jpype.JPackage('jline').lang.processes.Exp.fitMean(mean))

    fit_mean = fitMean



class Erlang(Markovian):
    """
    Erlang distribution for controlled-variance service times.
    
    An Erlang distribution is the sum of k independent exponential random variables,
    providing a way to model service times with coefficient of variation less than 1.
    
    Args:
        rate (float): Rate parameter for each phase.
        nphases (int): Number of phases (k).
        
    Examples:
        >>> Erlang(3.0, 2)  # 2-phase Erlang with rate 3.0 per phase
        >>> Erlang.fitMeanAndSCV(1.0, 0.5)  # Fit to mean=1.0, SCV=0.5
    """
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            self.obj = args[0]
        else:
            rate = args[0]
            nphases = args[1]
            self.obj = jpype.JPackage('jline').lang.processes.Erlang(rate, nphases)

    def fitMeanAndSCV(mean, scv):
        """
        Create an Erlang distribution with specified mean and squared coefficient of variation.
        
        Args:
            mean (float): Mean service time > 0.
            scv (float): Squared coefficient of variation, must be ≤ 1.0 for Erlang.
            
        Returns:
            Erlang: An Erlang distribution matching the given mean and SCV.
        """
        return Erlang(jpype.JPackage('jline').lang.processes.Erlang.fitMeanAndSCV(mean, scv))

    fit_mean_and_scv = fitMeanAndSCV

    def fitMeanAndOrder(mean, order):
        """
        Create an Erlang distribution with specified mean and number of phases.
        
        Args:
            mean (float): Mean service time > 0.
            order (int): Number of phases (k ≥ 1).
            
        Returns:
            Erlang: An Erlang distribution with k phases and the specified mean.
        """
        return Erlang(jpype.JPackage('jline').lang.processes.Erlang.fitMeanAndOrder(mean, order))

    fit_mean_and_order = fitMeanAndOrder


class Gamma(ContinuousDistribution):
    """
    Gamma distribution for flexible service time modeling.
    
    A continuous probability distribution that generalizes the exponential
    and Erlang distributions, allowing for flexible coefficient of variation.
    
    Args:
        shape (float): Shape parameter (α).
        scale (float): Scale parameter (β).
        
    Examples:
        >>> Gamma(2.0, 0.5)  # Gamma with shape=2.0, scale=0.5
        >>> Gamma.fitMeanAndSCV(1.0, 2.0)  # Fit to mean=1.0, SCV=2.0
    """
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            self.obj = args[0]
        else:
            shape = args[0]
            scale = args[1]
            self.obj = jpype.JPackage('jline').lang.processes.Gamma(shape, scale)

    def fitMeanAndSCV(mean, scv):
        return Gamma(jpype.JPackage('jline').lang.processes.Gamma.fitMeanAndSCV(mean, scv))

    fit_mean_and_scv = fitMeanAndSCV


class Geometric(DiscreteDistribution):
    """
    Geometric distribution for discrete waiting times.
    
    Represents the number of trials needed to get the first success
    in a sequence of independent Bernoulli trials.
    
    Args:
        prob (float): Success probability for each trial.
        
    Examples:
        >>> Geometric(0.1)  # Geometric with success probability 0.1
    """
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            if isinstance(args[0], jpype.JPackage('jline').lang.processes.DiscreteDistribution):
                self.obj = args[0]
            else:
                prob = args[0]
                self.obj = jpype.JPackage('jline').lang.processes.Geometric(prob)


class HyperExp(Markovian):
    """
    Hyper-exponential distribution for high-variability service times.
    
    A mixture of exponential distributions that can model service times
    with coefficient of variation greater than 1, commonly used when
    service times are highly variable.
    
    Args:
        p (float): Probability of selecting the first exponential component.
        lambda1 (float): Rate of the first exponential component.
        lambda2 (float): Rate of the second exponential component.
        
    Examples:
        >>> HyperExp(0.3, 1.0, 5.0)  # Mixture with prob=0.3, rates 1.0 and 5.0
        >>> HyperExp.fitMeanAndSCV(1.0, 4.0)  # Fit to mean=1.0, SCV=4.0
    """
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            if hasattr(args[0], '__module__') and 'jpype' in args[0].__module__:
                self.obj = args[0]
            elif hasattr(args[0], '__class__') and 'jline.lang.processes.HyperExp' in str(args[0].__class__):
                self.obj = args[0]
            else:
                self.obj = jpype.JPackage('jline').lang.processes.HyperExp(0.5, args[0])
        elif len(args) == 2:
            p = args[0]
            lambda_rate = args[1]
            self.obj = jpype.JPackage('jline').lang.processes.HyperExp(p, lambda_rate)
        else:
            p = args[0]
            lambda1 = args[1]
            lambda2 = args[2]
            self.obj = jpype.JPackage('jline').lang.processes.HyperExp(p, lambda1, lambda2)

    @staticmethod
    def fitMeanAndSCV(mean, scv):
        return HyperExp(jpype.JPackage('jline').lang.processes.HyperExp.fitMeanAndSCV(mean, scv))

    fit_mean_and_scv = fitMeanAndSCV

    @staticmethod
    def fitMeanAndSCVBalanced(mean, scv):
        """
        Create a balanced hyperexponential distribution with specified mean and SCV.
        
        A balanced hyperexponential uses equal probabilities for both phases,
        which can provide better numerical stability than the general case.
        
        Args:
            mean (float): Mean service time > 0.
            scv (float): Squared coefficient of variation, must be ≥ 1.0 for hyperexponential.
            
        Returns:
            HyperExp: A balanced hyperexponential distribution matching the given mean and SCV.
        """
        return HyperExp(jpype.JPackage('jline').lang.processes.HyperExp.fitMeanAndSCVBalanced(mean, scv))

    fit_mean_and_scv_balanced = fitMeanAndSCVBalanced


class Immediate(Distribution):
    def __init__(self):
        super().__init__()
        self.obj = jpype.JPackage('jline').lang.processes.Immediate()


class Lognormal(ContinuousDistribution):
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            self.obj = args[0]
        else:
            mu = args[0]
            sigma = args[1]
            self.obj = jpype.JPackage('jline').lang.processes.Lognormal(mu, sigma)

    def fitMeanAndSCV(mean, scv):
        return Lognormal(jpype.JPackage('jline').lang.processes.Lognormal.fitMeanAndSCV(mean, scv))

    fit_mean_and_scv = fitMeanAndSCV


class Normal(ContinuousDistribution):
    """
    Normal (Gaussian) distribution.

    Args:
        mu (float): Mean of the distribution.
        sigma (float): Standard deviation of the distribution.

    Examples:
        >>> normal = Normal(0.0, 1.0)  # Standard normal
        >>> normal.getMean()
        0.0
        >>> normal.getStd()
        1.0
    """
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            self.obj = args[0]
        else:
            mu = args[0]
            sigma = args[1]
            self.obj = jpype.JPackage('jline').lang.processes.Normal(mu, sigma)

    @staticmethod
    def fitMean(mean):
        """Create a Normal distribution with specified mean and unit variance."""
        return Normal(jpype.JPackage('jline').lang.processes.Normal.fitMean(mean))

    @staticmethod
    def fitMeanAndStd(mean, std):
        """Create a Normal distribution with specified mean and standard deviation."""
        return Normal(jpype.JPackage('jline').lang.processes.Normal.fitMeanAndStd(mean, std))

    @staticmethod
    def fitMeanAndVar(mean, var):
        """Create a Normal distribution with specified mean and variance."""
        return Normal(jpype.JPackage('jline').lang.processes.Normal.fitMeanAndVar(mean, var))

    def getMean(self):
        """Get the mean of the distribution."""
        return float(self.obj.getMean())

    def getStd(self):
        """Get the standard deviation of the distribution."""
        return float(self.obj.getStd())

    def getVar(self):
        """Get the variance of the distribution."""
        return float(self.obj.getVar())

    def getSCV(self):
        """Get the squared coefficient of variation."""
        return float(self.obj.getSCV())

    def getSkewness(self):
        """Get the skewness (always 0 for normal)."""
        return float(self.obj.getSkewness())

    def evalCDF(self, x):
        """Evaluate the cumulative distribution function at x."""
        return float(self.obj.evalCDF(x))

    def evalPDF(self, x):
        """Evaluate the probability density function at x."""
        return float(self.obj.evalPDF(x))

    def evalLST(self, s):
        """Evaluate the Laplace-Stieltjes transform at s."""
        return float(self.obj.evalLST(s))

    def sample(self, n, random=None):
        """
        Generate random samples from the distribution.

        Args:
            n (int): Number of samples to generate.
            random (Random, optional): Java Random instance for reproducibility.

        Returns:
            numpy.ndarray: Array of n samples.
        """
        import numpy as np
        if random is None:
            return np.array(self.obj.sample(n))
        else:
            return np.array(self.obj.sample(n, random))

    # Snake_case aliases
    fit_mean = fitMean
    fit_mean_and_std = fitMeanAndStd
    fit_mean_and_var = fitMeanAndVar
    get_mean = getMean
    get_std = getStd
    get_var = getVar
    get_scv = getSCV
    get_skewness = getSkewness
    eval_cdf = evalCDF
    eval_pdf = evalPDF
    eval_lst = evalLST


class MAP(Markovian):
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            self.obj = args[0]
        else:
            D0 = args[0]
            D1 = args[1]
            self.obj = jpype.JPackage('jline').lang.processes.MAP(jlineMatrixFromArray(D0), jlineMatrixFromArray(D1))

    def toPH(self):
        self.obj.toPH()

    @staticmethod
    def rand(order=None):
        """
        Create a random MAP distribution.

        Args:
            order (int, optional): The order (number of phases) of the MAP.
                                   Default is 2 if not specified.

        Returns:
            MAP: A randomly generated MAP distribution.
        """
        if order is None:
            return MAP(jpype.JPackage('jline').lang.processes.MAP.rand())
        else:
            return MAP(jpype.JPackage('jline').lang.processes.MAP.rand(int(order)))

    def setMean(self, new_mean):
        """
        Scale the MAP distribution to have a specified mean.

        Args:
            new_mean (float): The desired mean inter-arrival/service time.

        Returns:
            MAP: Self, to allow method chaining.
        """
        self.obj.setMean(float(new_mean))
        return self

    set_mean = setMean

    def getMean(self):
        """
        Get the mean of the MAP distribution.

        Returns:
            float: The mean inter-arrival/service time.
        """
        return float(self.obj.getMean())

    get_mean = getMean

    def toMMDP(self):
        """
        Convert MAP to MMDP (deterministic representation for fluid queues).

        Converts this Markovian Arrival Process to a Markov-Modulated
        Deterministic Process suitable for fluid queue analysis.

        Returns:
            MMDP: MMDP representation of this MAP
        """
        return MMDP.fromMAP(self)

    to_mmdp = toMMDP


class MarkedMAP(Markovian):
    """
    Marked Markovian Arrival Process (MarkedMAP/MMAP).

    Represents a Markovian arrival process where arrivals are marked/tagged
    with different types. Uses the M3A representation format.

    Args:
        D: MatrixCell or list containing {D0, D1, D11, D12, ..., D1K}
        K: Number of marking types (optional if using object constructor)
    """

    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            # Existing Java object
            self.obj = args[0]
        elif len(args) == 2:
            # D matrices and K
            D = args[0]
            K = args[1]
            # Convert Python list/arrays to Java MatrixCell
            from jline.util.matrix import MatrixCell
            mc = MatrixCell()
            for i, mat in enumerate(D):
                mc.set(i, jlineMatrixFromArray(mat))
            self.obj = jpype.JPackage('jline').lang.processes.MarkedMAP(mc)
        else:
            raise ValueError("MarkedMAP requires either 1 or 2 arguments")

    def getNumberOfTypes(self):
        """
        Get the number of marking types.

        Returns:
            int: Number of distinct mark types
        """
        return int(self.obj.getNumberOfTypes())

    get_number_of_types = getNumberOfTypes

    def normalize(self):
        """
        Normalize the MarkedMAP so that D0+sum(Di) forms a proper generator.

        Returns:
            MarkedMAP: Self, to allow method chaining
        """
        self.obj.normalize()
        return self

    @staticmethod
    def rand(order=2, num_types=2):
        """
        Create a random MarkedMAP distribution.

        Args:
            order (int): Number of phases. Default is 2.
            num_types (int): Number of marking types. Default is 2.

        Returns:
            MarkedMAP: A randomly generated MarkedMAP distribution
        """
        return MarkedMAP(jpype.JPackage('jline').lang.processes.MarkedMAP.rand(int(order), int(num_types)))


class BMAP(MarkedMAP):
    """
    Batch Markovian Arrival Process (BMAP).

    A point process where arrivals occur in batches. Uses standard BMAP representation:
    - D0: infinitesimal generator for transitions without arrivals
    - D1: rate matrix for transitions generating 1 arrival
    - D2: rate matrix for transitions generating 2 arrivals
    - Dk: rate matrix for transitions generating k arrivals

    BMAP extends MarkedMAP where each mark k represents a batch size k.

    Args:
        Can be initialized in several ways:
        - BMAP(java_obj): Wrap existing Java BMAP object
        - BMAP(D0, D1, D2, ...): Create from D matrices as separate numpy arrays
        - BMAP([D0, D1, D2, ...]): Create from list of D matrices
    """

    def __init__(self, *args):
        Markovian.__init__(self)  # Skip MarkedMAP init, go straight to Markovian

        if len(args) == 1:
            if isinstance(args[0], jpype.JClass('jline.lang.processes.BMAP')):
                # Existing Java BMAP object
                self.obj = args[0]
            elif isinstance(args[0], list):
                # List of D matrices [D0, D1, D2, ...]
                D_list = args[0]
                # Convert to Java varargs format
                D0 = jlineMatrixFromArray(D_list[0])
                Dk = [jlineMatrixFromArray(d) for d in D_list[1:]]
                # Call Java constructor with varargs
                self.obj = jpype.JPackage('jline').lang.processes.BMAP(D0, *Dk)
            else:
                raise ValueError("Single argument must be Java BMAP object or list of matrices")
        elif len(args) >= 2:
            # Multiple D matrices as separate arguments: BMAP(D0, D1, D2, ...)
            D0 = jlineMatrixFromArray(args[0])
            Dk = [jlineMatrixFromArray(d) for d in args[1:]]
            self.obj = jpype.JPackage('jline').lang.processes.BMAP(D0, *Dk)
        else:
            raise ValueError("BMAP requires at least 1 argument")

    @staticmethod
    def fromMAPWithBatchPMF(D0, D1, batch_sizes, pmf):
        """
        Create BMAP from a base MAP and batch size distribution.

        Given a base MAP (D0, D1) for inter-batch arrivals and a batch size
        distribution, constructs the BMAP by scaling: Dk = D1 * pmf[k-1]

        Args:
            D0 (array): Base MAP's D0 matrix (transitions without batch arrivals)
            D1 (array): Base MAP's D1 matrix (inter-batch arrival transitions)
            batch_sizes (list): Array of batch sizes (e.g., [1, 2, 4, 8])
            pmf (list): Probability mass function for batch sizes

        Returns:
            BMAP: Constructed BMAP from the base MAP and batch distribution
        """
        D0_java = jlineMatrixFromArray(D0)
        D1_java = jlineMatrixFromArray(D1)
        batch_sizes_java = jpype.JArray(jpype.JInt)(batch_sizes)
        pmf_java = jpype.JArray(jpype.JDouble)(pmf)

        java_bmap = jpype.JPackage('jline').lang.processes.BMAP.fromMAPWithBatchPMF(
            D0_java, D1_java, batch_sizes_java, pmf_java)
        return BMAP(java_bmap)

    from_map_with_batch_pmf = fromMAPWithBatchPMF

    def getMaxBatchSize(self):
        """
        Get the maximum batch size supported by this BMAP.

        Returns:
            int: Maximum batch size k where Dk is defined
        """
        return int(self.obj.getMaxBatchSize())

    get_max_batch_size = getMaxBatchSize

    def getMeanBatchSize(self):
        """
        Get the mean batch size (expected number of arrivals per batch).

        Computed as weighted average of batch sizes by their rates:
        E[batch size] = sum(k * rate_k) / sum(rate_k)

        Returns:
            float: Mean batch size
        """
        return float(self.obj.getMeanBatchSize())

    get_mean_batch_size = getMeanBatchSize

    def getBatchRates(self):
        """
        Get the arrival rates for each batch size.

        Returns:
            ndarray: Array where result[k-1] is the rate of batch size k arrivals
        """
        rates_java = self.obj.getBatchRates()
        return np.array([rates_java[i] for i in range(len(rates_java))])

    get_batch_rates = getBatchRates

    def getBatchMatrix(self, batch_size):
        """
        Get the D_k matrix for batch size k.

        Args:
            batch_size (int): The batch size (1, 2, 3, ...)

        Returns:
            ndarray: D_k matrix for the specified batch size, or None if out of range
        """
        mat = self.obj.getBatchMatrix(int(batch_size))
        if mat is None:
            return None
        return jlineMatrixToArray(mat)

    get_batch_matrix = getBatchMatrix

    @staticmethod
    def rand(order=2, max_batch_size=3):
        """
        Create a random BMAP distribution.

        Args:
            order (int): Number of phases. Default is 2.
            max_batch_size (int): Maximum batch size. Default is 3.

        Returns:
            BMAP: A randomly generated BMAP distribution
        """
        return BMAP(jpype.JPackage('jline').lang.processes.BMAP.rand(int(order), int(max_batch_size)))


class ME(Markovian):
    """
    Matrix Exponential (ME) distribution.

    The ME distribution is a generalization of Phase-Type (PH) distributions where
    the initial probability vector alpha can have entries outside [0,1] and the
    matrix parameter A has arbitrary structure (not necessarily a valid sub-generator).

    ME distributions are characterized by:
    - alpha: Initial vector (may have negative entries or sum != 1)
    - A: Matrix parameter (must have all eigenvalues with negative real parts)

    Args:
        alpha: Initial probability vector as Python list or numpy array
        A: Matrix parameter as 2D Python list or numpy array

    Examples:
        >>> # Create ME distribution with 2 phases
        >>> alpha = [0.3, 0.7]
        >>> A = [[-2.0, 1.5], [0.5, -3.0]]
        >>> me = ME(alpha, A)
        >>> me.getMean()

        >>> # Create ME from exponential distribution
        >>> me_exp = ME.fromExp(2.0)  # rate = 2.0

        >>> # Create ME from Erlang distribution
        >>> me_erlang = ME.fromErlang(3, 1.0)  # k=3, rate=1.0

        >>> # Fit ME to moments
        >>> moments = [1.0, 2.0, 6.0]  # Exponential moments
        >>> me_fitted = ME.fitMoments(moments)
    """

    def __init__(self, *args):
        super().__init__()

        if len(args) == 1:
            # Existing Java ME object
            self.obj = args[0]
        elif len(args) == 2:
            # Create from alpha and A matrices
            alpha = args[0]
            A = args[1]

            # Convert Python lists/arrays to Java Matrix
            alpha_matrix = jlineMatrixFromArray(alpha)
            A_matrix = jlineMatrixFromArray(A)

            # Create Java ME object
            self.obj = jpype.JPackage('jline').lang.processes.ME(alpha_matrix, A_matrix)
        else:
            raise ValueError("ME constructor requires either 1 argument (Java object) or 2 arguments (alpha, A)")

    def getAlpha(self):
        """
        Get the initial probability vector alpha.

        Returns:
            numpy.ndarray: The alpha vector as a 1D numpy array
        """
        result = jlineMatrixToArray(self.obj.getAlpha())
        return result.flatten() if result is not None else None

    def getA(self):
        """
        Get the matrix parameter A.

        Returns:
            numpy.ndarray: The A matrix as a 2D numpy array
        """
        return jlineMatrixToArray(self.obj.getA())

    def getNumberOfPhases(self):
        """
        Get the number of phases in the ME representation.

        Returns:
            int: Number of phases
        """
        return int(self.obj.getNumberOfPhases())

    get_alpha = getAlpha
    get_a = getA
    get_number_of_phases = getNumberOfPhases

    @staticmethod
    def fromExp(rate):
        """
        Create ME distribution from exponential distribution.

        Args:
            rate (float): Rate parameter (lambda)

        Returns:
            ME: ME distribution equivalent to Exp(rate)
        """
        return ME(jpype.JPackage('jline').lang.processes.ME.fromExp(float(rate)))

    @staticmethod
    def fromErlang(k, rate):
        """
        Create ME distribution from Erlang distribution.

        Args:
            k (int): Number of phases
            rate (float): Rate parameter for each phase

        Returns:
            ME: ME distribution equivalent to Erlang(k, rate)
        """
        return ME(jpype.JPackage('jline').lang.processes.ME.fromErlang(int(k), float(rate)))

    @staticmethod
    def fromHyperExp(p, rates):
        """
        Create ME distribution from HyperExponential distribution.

        Args:
            p: Array of probabilities for each branch (Python list or numpy array)
            rates: Array of rates for each branch (Python list or numpy array)

        Returns:
            ME: ME distribution equivalent to HyperExp(p, rates)
        """
        # Convert to Java double[] arrays
        p_list = list(p) if isinstance(p, np.ndarray) else p
        rates_list = list(rates) if isinstance(rates, np.ndarray) else rates
        p_java = jpype.JArray(jpype.JDouble)(p_list)
        rates_java = jpype.JArray(jpype.JDouble)(rates_list)
        return ME(jpype.JPackage('jline').lang.processes.ME.fromHyperExp(p_java, rates_java))

    @staticmethod
    def fitMoments(moments):
        """
        Create ME distribution by fitting the given moments.
        Uses BuTools MEFromMoments algorithm.

        Args:
            moments: Array of moments (Python list or numpy array)
                    Requires 2*M-1 moments for order M ME

        Returns:
            ME: ME distribution matching the given moments
        """
        # Convert to Java double[] array
        moments_list = list(moments) if isinstance(moments, np.ndarray) else moments
        moments_java = jpype.JArray(jpype.JDouble)(moments_list)
        return ME(jpype.JPackage('jline').lang.processes.ME.fitMoments(moments_java))

    from_exp = fromExp
    from_erlang = fromErlang
    from_hyper_exp = fromHyperExp
    fit_moments = fitMoments


class RAP(Markovian):
    """
    Rational Arrival Process (RAP) distribution.

    RAP is a generalization of the Markovian Arrival Process (MAP) where the matrices
    H0 and H1 represent hidden and visible transitions respectively, but with relaxed
    constraints compared to MAP.

    RAP distributions are characterized by:
    - H0: matrix of hidden transition rates (transitions without arrivals)
    - H1: matrix of visible transition rates (transitions with arrivals)
    - H0 + H1 must form a valid infinitesimal generator (row sums = 0)
    - All eigenvalues of H0 must have negative real parts
    - Dominant eigenvalue of H0 must be negative and real

    The marginal distribution of inter-arrival times is a Matrix Exponential (ME) distribution.

    Args:
        H0: Hidden transition matrix as 2D Python list or numpy array
        H1: Visible transition matrix as 2D Python list or numpy array

    Examples:
        >>> # Create RAP distribution with 2 phases
        >>> H0 = [[-2.0, 1.0], [0.5, -1.5]]
        >>> H1 = [[0.5, 0.5], [0.5, 0.5]]
        >>> rap = RAP(H0, H1)
        >>> rap.getMean()

        >>> # Create RAP from Poisson process
        >>> rap_poisson = RAP.fromPoisson(2.0)  # rate = 2.0

        >>> # Create RAP from Erlang renewal process
        >>> rap_erlang = RAP.fromErlang(2, 1.0)  # k=2, rate=1.0

        >>> # Convert MAP to RAP
        >>> map_dist = MAP(D0, D1)
        >>> rap_from_map = RAP.fromMAP(map_dist)
    """

    def __init__(self, *args):
        super().__init__()

        if len(args) == 1:
            # Existing Java RAP object or MAP object
            if isinstance(args[0], MAP):
                # Convert MAP to RAP
                map_obj = args[0]
                self.obj = jpype.JPackage('jline').lang.processes.RAP.fromMAP(map_obj.obj)
            else:
                # Existing Java RAP object
                self.obj = args[0]
        elif len(args) == 2:
            # Create from H0 and H1 matrices
            H0 = args[0]
            H1 = args[1]

            # Convert Python lists/arrays to Java Matrix
            H0_matrix = jlineMatrixFromArray(H0)
            H1_matrix = jlineMatrixFromArray(H1)

            # Create Java RAP object
            self.obj = jpype.JPackage('jline').lang.processes.RAP(H0_matrix, H1_matrix)
        else:
            raise ValueError("RAP constructor requires either 1 argument (Java object/MAP) or 2 arguments (H0, H1)")

    def getH0(self):
        """
        Get the H0 matrix (hidden transition rates).

        Returns:
            numpy.ndarray: The H0 matrix as a 2D numpy array
        """
        return jlineMatrixToArray(self.obj.getH0())

    def getH1(self):
        """
        Get the H1 matrix (visible transition rates).

        Returns:
            numpy.ndarray: The H1 matrix as a 2D numpy array
        """
        return jlineMatrixToArray(self.obj.getH1())

    def getNumberOfPhases(self):
        """
        Get the number of phases in the RAP representation.

        Returns:
            int: Number of phases
        """
        return int(self.obj.getNumberOfPhases())

    def getRate(self):
        """
        Get the arrival rate (lambda) of the RAP.

        Returns:
            float: Arrival rate
        """
        return float(self.obj.getRate())

    get_h0 = getH0
    get_h1 = getH1
    get_number_of_phases = getNumberOfPhases
    get_rate = getRate

    @staticmethod
    def fromPoisson(rate):
        """
        Create RAP from exponential renewal process (Poisson process).

        Args:
            rate (float): Arrival rate (lambda)

        Returns:
            RAP: RAP distribution representing a Poisson process
        """
        return RAP(jpype.JPackage('jline').lang.processes.RAP.fromPoisson(float(rate)))

    @staticmethod
    def fromErlang(k, rate):
        """
        Create RAP from Erlang renewal process.

        Args:
            k (int): Number of phases
            rate (float): Rate parameter for each phase

        Returns:
            RAP: RAP distribution representing an Erlang renewal process
        """
        return RAP(jpype.JPackage('jline').lang.processes.RAP.fromErlang(int(k), float(rate)))

    @staticmethod
    def fromMAP(map_dist):
        """
        Create RAP from a Markovian Arrival Process (MAP).
        This shows that MAP is a special case of RAP.

        Args:
            map_dist (MAP): MAP distribution instance

        Returns:
            RAP: RAP distribution equivalent to the given MAP
        """
        if isinstance(map_dist, MAP):
            return RAP(jpype.JPackage('jline').lang.processes.RAP.fromMAP(map_dist.obj))
        else:
            raise ValueError("Input must be a MAP object")

    from_poisson = fromPoisson
    from_erlang = fromErlang
    from_map = fromMAP


class MMPP2(Markovian):
    """
    Two-state Markov Modulated Poisson Process (MMPP/2).

    Represents a Poisson process where the arrival rate is modulated
    by a two-state Markov chain, allowing for correlated arrivals.

    Args:
        lambda0: Arrival rate in state 0
        lambda1: Arrival rate in state 1
        sigma0: Transition rate from state 0 to state 1
        sigma1: Transition rate from state 1 to state 0
    """
    def __init__(self, lambda0, lambda1, sigma0, sigma1):
        super().__init__()
        self.obj = jpype.JPackage('jline').lang.processes.MMPP2(lambda0, lambda1, sigma0, sigma1)

    @staticmethod
    def fitCentral(mean, var, skew):
        """
        Create an MMPP2 distribution matching the first three central moments.

        Args:
            mean (float): Mean (first moment).
            var (float): Variance (second central moment).
            skew (float): Skewness (third standardized moment).

        Returns:
            MMPP2: An MMPP2 distribution matching the given moments.
        """
        return MMPP2(jpype.JPackage('jline').lang.processes.MMPP2.fitCentral(mean, var, skew))

    fit_central = fitCentral

    def toMMDP(self):
        """
        Convert MMPP2 to MMDP2 (deterministic representation for fluid queues).

        Converts this Markov-Modulated Poisson Process to a Markov-Modulated
        Deterministic Process suitable for fluid queue analysis. Uses the same
        parameterization where Poisson rates become deterministic flow rates.

        Returns:
            MMDP2: MMDP2 representation with the same parameters.
        """
        return MMDP2(self.obj.toMMDP())

    to_mmdp = toMMDP


class MMDP(Markovian):
    """
    Markov-Modulated Deterministic Process (MMDP) for fluid queue modeling.

    Models fluid flow with deterministic rates modulated by a background
    Markov chain. Suitable for arrival and service processes in Markovian
    fluid queues analyzed by the mfq method of SolverFLD.

    The (Q, R) parameterization follows BUTools conventions:
    - Q: Generator matrix of the modulating CTMC (row sums = 0)
    - R: Diagonal matrix of deterministic rates per state

    MMDP is the deterministic analogue of MMPP:
    - MMPP: Poisson arrival rates modulated by a Markov chain
    - MMDP: Deterministic rates modulated by a Markov chain

    Args:
        Q: n×n generator matrix (row sums = 0), or existing Java MMDP object
        R: n×n diagonal matrix of rates, or n-vector (optional if Q is Java object)

    Examples:
        >>> # Create 2-state MMDP with different flow rates
        >>> Q = [[-0.5, 0.5], [0.3, -0.3]]
        >>> R = [[2.0, 0], [0, 5.0]]  # or just [2.0, 5.0]
        >>> mmdp = MMDP(Q, R)
        >>> mean_rate = mmdp.getMeanRate()  # Stationary mean flow rate
    """

    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            # Existing Java MMDP object
            self.obj = args[0]
        elif len(args) == 2:
            Q = args[0]
            R = args[1]
            self.obj = jpype.JPackage('jline').lang.processes.MMDP(
                jlineMatrixFromArray(Q), jlineMatrixFromArray(R))
        else:
            raise ValueError("MMDP requires either 1 argument (Java object) or 2 arguments (Q, R)")

    def Q(self):
        """
        Get the generator matrix Q.

        Returns:
            numpy.ndarray: n×n generator matrix of the modulating CTMC
        """
        return jlineMatrixToArray(self.obj.Q())

    def R(self):
        """
        Get the rate matrix R (diagonal).

        Returns:
            numpy.ndarray: n×n diagonal matrix of deterministic rates
        """
        return jlineMatrixToArray(self.obj.R())

    def r(self):
        """
        Get the rate vector (diagonal of R).

        Returns:
            numpy.ndarray: n-vector of deterministic rates per phase
        """
        return jlineMatrixToArray(self.obj.r()).flatten()

    def getNumberOfPhases(self):
        """
        Get the number of phases in the modulating CTMC.

        Returns:
            int: Number of phases
        """
        return int(self.obj.getNumberOfPhases())

    def getMeanRate(self):
        """
        Compute the stationary mean rate.

        For MMDP: E[r] = π * diag(R), where π is the stationary distribution
        of the modulating CTMC with generator Q.

        Returns:
            float: Stationary mean deterministic rate
        """
        return float(self.obj.getMeanRate())

    def getSCV(self):
        """
        Compute the squared coefficient of variation of rates.

        Computes Var[r]/E[r]^2 where expectation is over the
        stationary distribution of the modulating CTMC.

        Returns:
            float: Squared coefficient of variation
        """
        return float(self.obj.getSCV())

    def getProcess(self):
        """
        Get the (Q, R) representation as a tuple.

        Returns:
            tuple: (Q, R) matrices
        """
        return (self.Q(), self.R())

    get_number_of_phases = getNumberOfPhases
    get_mean_rate = getMeanRate
    get_scv = getSCV
    get_process = getProcess

    @staticmethod
    def fromMAP(map_dist):
        """
        Convert a MAP to MMDP (deterministic representation).

        Converts a Markovian Arrival Process to a Markov-Modulated
        Deterministic Process by extracting the full generator and
        using row sums of D1 as the deterministic rates.

        Args:
            map_dist (MAP): MAP distribution to convert

        Returns:
            MMDP: MMDP representation of the MAP
        """
        if isinstance(map_dist, MAP):
            return MMDP(jpype.JPackage('jline').lang.processes.MMDP.fromMAP(map_dist.obj))
        else:
            raise ValueError("Input must be a MAP object")

    @staticmethod
    def fromMMPP2(lambda0, lambda1, sigma0, sigma1):
        """
        Create MMDP from MMPP2 parameters.

        Creates a 2-state MMDP using the same parameterization as MMPP2.

        Args:
            lambda0: Rate in state 0
            lambda1: Rate in state 1
            sigma0: Transition rate from state 0 to state 1
            sigma1: Transition rate from state 1 to state 0

        Returns:
            MMDP: 2-state MMDP
        """
        return MMDP(jpype.JPackage('jline').lang.processes.MMDP.fromMMPP2(
            float(lambda0), float(lambda1), float(sigma0), float(sigma1)))

    @staticmethod
    def isFeasible(Q, R):
        """
        Check if (Q, R) defines a valid MMDP.

        Requirements:
        - Q must be a valid generator (square, row sums = 0, proper signs)
        - R must be diagonal with non-negative entries
        - Q and R must have compatible dimensions

        Args:
            Q: Generator matrix
            R: Rate matrix

        Returns:
            bool: True if valid MMDP, False otherwise
        """
        return bool(jpype.JPackage('jline').lang.processes.MMDP.isFeasible(
            jlineMatrixFromArray(Q), jlineMatrixFromArray(R)))

    from_map = fromMAP
    from_mmpp2 = fromMMPP2
    is_feasible = isFeasible


class MMDP2(MMDP):
    """
    2-state Markov-Modulated Deterministic Process (MMDP/2).

    A specialized MMDP with exactly 2 phases, using a convenient
    parameterization analogous to MMPP2.

    Parameterization:
        r0, r1: Deterministic rates in states 0 and 1
        sigma0: Transition rate from state 0 to state 1
        sigma1: Transition rate from state 1 to state 0

    The generator matrix is:
        Q = [[-sigma0, sigma0], [sigma1, -sigma1]]

    The rate matrix is:
        R = diag([r0, r1])

    Args:
        r0: Deterministic rate in state 0, or existing Java MMDP object
        r1: Deterministic rate in state 1 (optional if r0 is Java object)
        sigma0: Transition rate from state 0 to state 1
        sigma1: Transition rate from state 1 to state 0

    Examples:
        >>> mmdp2 = MMDP2(2.0, 5.0, 0.5, 0.3)
        >>> mean_rate = mmdp2.getMeanRate()  # (2.0*0.3 + 5.0*0.5) / (0.5 + 0.3)
    """

    def __init__(self, *args):
        Markovian.__init__(self)  # Skip MMDP.__init__
        if len(args) == 1:
            # Existing Java MMDP2 object
            self.obj = args[0]
        elif len(args) == 4:
            r0 = args[0]
            r1 = args[1]
            sigma0 = args[2]
            sigma1 = args[3]
            # Build Q and R matrices for 2-state case
            Q = [[-sigma0, sigma0], [sigma1, -sigma1]]
            R = [[r0, 0], [0, r1]]
            self.obj = jpype.JPackage('jline').lang.processes.MMDP(
                jlineMatrixFromArray(Q), jlineMatrixFromArray(R))
        else:
            raise ValueError("MMDP2 requires 1 argument (Java object) or 4 arguments (r0, r1, sigma0, sigma1)")

    def getNumberOfPhases(self):
        """
        Get the number of phases.

        Returns:
            int: Always 2 for MMDP2
        """
        return 2

    get_number_of_phases = getNumberOfPhases


class PH(Markovian):
    """
    Phase-type (PH) distribution.
    
    Represents a phase-type distribution defined by initial probability
    vector and sub-generator matrix. PH distributions model the time
    until absorption in a finite Markov chain.
    
    Args:
        alpha: Initial probability vector or existing PH object
        subgen: Sub-generator matrix (if alpha is vector)
    """
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            self.obj = args[0]
        else:
            alpha = args[0]
            subgen = args[1]
            self.obj = jpype.JPackage('jline').lang.processes.PH(jlineMatrixFromArray(alpha),
                                                                     jlineMatrixFromArray(subgen))

    @staticmethod
    def fitCentral(mean, var, skew):
        """
        Create a PH distribution matching the first three central moments.

        Args:
            mean (float): Mean (first moment).
            var (float): Variance (second central moment).
            skew (float): Skewness (third standardized moment).

        Returns:
            PH: A PH distribution matching the given moments.
        """
        return PH(jpype.JPackage('jline').lang.processes.PH.fitCentral(mean, var, skew))

    fit_central = fitCentral


class Pareto(ContinuousDistribution):
    """
    Pareto distribution (power law distribution).
    
    Represents a Pareto distribution with specified shape and scale
    parameters, commonly used for modeling heavy-tailed phenomena.
    
    Args:
        shape: Shape parameter (alpha)
        scale: Scale parameter (x_min) or existing Pareto object
    """
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            self.obj = args[0]
        else:
            shape = args[0]
            scale = args[1]
            self.obj = jpype.JPackage('jline').lang.processes.Pareto(shape, scale)

    @staticmethod
    def fitMeanAndSCV(mean, scv):
        return Pareto(jpype.JPackage('jline').lang.processes.Pareto.fitMeanAndSCV(mean, scv))

    fit_mean_and_scv = fitMeanAndSCV


class Poisson(DiscreteDistribution):
    """
    Poisson distribution for count data.
    
    Represents a discrete probability distribution expressing the
    probability of a given number of events occurring in a fixed
    interval when events occur with a known constant rate.
    
    Args:
        rate: Expected number of events (lambda parameter)
    """
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            self.obj = args[0]
        else:
            rate = args[0]
            self.obj = jpype.JPackage('jline').lang.processes.Poisson(rate)


class Replayer(Distribution):
    """
    Empirical distribution that replays trace data from files.

    Replayer reads real-world trace data (time series) from files and reproduces
    the exact sequence of observations. This is essential for trace-driven
    simulation, workload characterization, and model validation against real data.

    Supports lazy loading of trace data to manage memory efficiently. Provides
    statistical analysis of traces and distribution fitting.

    Args:
        data: Either a filename (str) or a Distribution object to wrap

    Example:
        >>> replayer = Replayer('arrival_times.txt')
        >>> mean = replayer.getMean()
        >>> replayer.load()  # Load into memory
        >>> sample = replayer.sample()  # Get next value
        >>> replayer.unload()  # Free memory
        >>> exp_dist = replayer.fitExp()  # Fit exponential
        >>> cox_dist = replayer.fitCoxian()  # Fit Cox2
    """

    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            if isinstance(args[0], Distribution):
                self.obj = args[0]
                self._data = None
                self._cursample = 0
                self._filename = None
            else:
                filename = args[0]
                self.obj = jpype.JPackage('jline').lang.processes.Replayer(filename)
                self._data = None  # Lazy loading
                self._cursample = 0
                self._filename = filename

    def load(self):
        """
        Load trace data from file into memory.

        This reads all values from the trace file (one value per line) and
        stores them in memory for efficient sampling. After loading, getters
        like getMean(), getSCV(), etc. can access the data without re-reading.

        Raises:
            ValueError: If the Replayer was not constructed from a file
            FileNotFoundError: If the trace file cannot be found
            ValueError: If the file contains non-numeric values

        Example:
            >>> replayer = Replayer('trace.txt')
            >>> replayer.load()  # Load all data
            >>> mean = replayer.getMean()
        """
        if self._filename is None:
            raise ValueError("Cannot load: Replayer not constructed from file")

        try:
            with open(self._filename, 'r') as f:
                self._data = []
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if line:  # Skip empty lines
                        try:
                            self._data.append(float(line))
                        except ValueError:
                            raise ValueError(
                                f"Invalid number in trace file '{self._filename}' "
                                f"at line {line_num}: {line}"
                            )
        except FileNotFoundError:
            raise FileNotFoundError(f"Trace file not found: {self._filename}")
        self._cursample = 0

    def unload(self):
        """
        Unload trace data from memory to free resources.

        After calling unload(), the trace data is discarded from memory.
        Subsequent calls to sample() or methods that need the data will
        trigger load() automatically.

        Example:
            >>> replayer = Replayer('trace.txt')
            >>> replayer.load()
            >>> # Use replayer...
            >>> replayer.unload()  # Free memory
        """
        self._data = None
        self._cursample = 0

    def getFileName(self):
        """
        Get the trace file name.

        Returns:
            str: Path to the trace file, or None if constructed from array
        """
        return self.obj.getFileName()

    def fitExp(self):
        """
        Fit exponential distribution to trace statistics.

        Creates an exponential distribution with mean matching the trace mean.
        This is useful for replacing empirical traces with analytical models.

        Returns:
            Exp: Exponential distribution with trace mean

        Example:
            >>> replayer = Replayer('trace.txt')
            >>> exp_dist = replayer.fitExp()  # Single parameter Exp(lambda)
        """
        return Exp.fitMean(self.getMean())

    def fitCoxian(self):
        """
        Fit Coxian-2 (Cox2) distribution to trace moments.

        Creates a Coxian-2 distribution matching the first three moments
        (mean, variance, skewness) of the trace. This provides a more
        flexible fit than exponential for non-exponential traces.

        Returns:
            Cox2: Coxian-2 distribution matching trace moments

        Example:
            >>> replayer = Replayer('trace.txt')
            >>> cox_dist = replayer.fitCoxian()  # 3-parameter Coxian
        """
        mean = self.getMean()
        scv = self.getSCV()
        skew = self.obj.getSkewness()  # Use Java method directly
        var = scv * mean * mean
        return Cox2.fitCentral(mean, var, skew)

    def sample(self):
        """
        Sample next value from trace (sequential replay).

        Returns successive trace values in order. When the end of trace is
        reached, raises an error. Use load() to explicitly load data, or
        this method will load automatically on first call.

        Returns:
            float: Next trace value

        Raises:
            ValueError: If end of trace data reached

        Example:
            >>> replayer = Replayer('trace.txt')
            >>> for i in range(10):
            ...     value = replayer.sample()  # Get next 10 values in order
        """
        if self._data is None:
            self.load()

        if self._cursample >= len(self._data):
            raise ValueError("End of trace data reached")

        value = self._data[self._cursample]
        self._cursample += 1
        return value

    def fitAPH(self):
        """
        Fit acyclic phase-type (APH) distribution to trace statistics.

        Uses the Java backend to fit an APH distribution to the trace.

        Returns:
            APH: Acyclic phase-type distribution
        """
        return APH(self.obj.fitAPH())

    # Snake_case aliases
    get_file_name = getFileName
    fit_exp = fitExp
    fit_coxian = fitCoxian
    fit_aph = fitAPH


class Uniform(ContinuousDistribution):
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            self.obj = args[0]
        else:
            minVal = args[0]
            maxVal = args[1]
            self.obj = jpype.JPackage('jline').lang.processes.Uniform(minVal, maxVal)


class Weibull(ContinuousDistribution):
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            self.obj = args[0]
        else:
            shape = args[0]
            scale = args[1]
            self.obj = jpype.JPackage('jline').lang.processes.Weibull(shape, scale)

    def fitMeanAndSCV(mean, scv):
        return Weibull(jpype.JPackage('jline').lang.processes.Weibull.fitMeanAndSCV(mean, scv))

    fit_mean_and_scv = fitMeanAndSCV


class Zipf(DiscreteDistribution):
    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            self.obj = args[0]
        else:
            s = args[0]
            n = args[1]
            self.obj = jpype.JPackage('jline').lang.processes.Zipf(s, n)


class Prior(Distribution):
    """
    Prior distribution for Bayesian-style parameter uncertainty modeling.

    Prior represents parameter uncertainty by specifying a discrete set of
    alternative distributions with associated probabilities. When used with
    setService or setArrival, it causes the Posterior solver to expand the
    model into a family of networks, one for each alternative.

    This is NOT a mixture distribution - each alternative represents a
    separate model realization with its associated prior probability.

    Args:
        distributions (list): List of Distribution objects representing alternatives.
        probabilities (list): List of probabilities (must sum to 1, be non-negative).

    Examples:
        >>> # Service time with uncertain rate
        >>> prior = Prior([Exp(1.0), Exp(2.0), Erlang(1.5, 2)], [0.4, 0.35, 0.25])
        >>> queue.setService(job_class, prior)
        >>>
        >>> # Solve with Posterior wrapper
        >>> from line_solver import SolverPosterior, SolverMVA
        >>> post = SolverPosterior(model, lambda m: SolverMVA(m))
        >>> avg_table = post.getAvgTable()  # Prior-weighted expectations
    """

    def __init__(self, *args):
        super().__init__()
        if len(args) == 1:
            # Existing Java object
            self.obj = args[0]
        elif len(args) == 2:
            distributions = args[0]
            probabilities = args[1]

            # Validate inputs
            if not isinstance(distributions, list):
                raise ValueError("distributions must be a list")
            if not isinstance(probabilities, list):
                probabilities = list(probabilities)
            if len(distributions) != len(probabilities):
                raise ValueError("Number of distributions must match number of probabilities")

            # Convert Python Distribution objects to Java
            java_list = jpype.JPackage('java').util.ArrayList()
            for dist in distributions:
                if hasattr(dist, 'obj'):
                    java_list.add(dist.obj)
                else:
                    java_list.add(dist)

            # Convert probabilities to Java double array
            prob_array = jpype.JArray(jpype.JDouble)(probabilities)

            # Create Java Prior object
            self.obj = jpype.JPackage('jline').lang.processes.Prior(java_list, prob_array)
        else:
            raise ValueError("Prior requires either 1 argument (Java object) or 2 arguments (distributions, probabilities)")

    def getNumAlternatives(self):
        """
        Get the number of alternative distributions.

        Returns:
            int: Number of alternatives in the prior.
        """
        return int(self.obj.getNumAlternatives())

    def getAlternative(self, idx):
        """
        Get the distribution at the specified index.

        Args:
            idx (int): 0-based index of the alternative.

        Returns:
            Distribution: The distribution at the given index.
        """
        java_dist = self.obj.getAlternative(idx)
        # Wrap in appropriate Python Distribution class
        return Distribution._wrap(java_dist)

    def getProbability(self, idx):
        """
        Get the probability of the alternative at the specified index.

        Args:
            idx (int): 0-based index of the alternative.

        Returns:
            float: The probability of the alternative.
        """
        return float(self.obj.getProbability(idx))

    def getProbabilities(self):
        """
        Get all probabilities as a list.

        Returns:
            list: List of all probabilities.
        """
        java_probs = self.obj.getProbabilities()
        return [float(p) for p in java_probs]

    def isPrior(self):
        """
        Check if this is a Prior distribution.

        Returns:
            bool: Always True for Prior instances.
        """
        return True

    @staticmethod
    def isPriorDistribution(dist):
        """
        Check if a distribution is a Prior.

        Args:
            dist: Distribution object to check.

        Returns:
            bool: True if dist is a Prior.
        """
        if isinstance(dist, Prior):
            return True
        if hasattr(dist, 'obj'):
            return jpype.JPackage('jline').lang.processes.Prior.isPriorDistribution(dist.obj)
        return False

    # Snake_case aliases
    get_num_alternatives = getNumAlternatives
    get_alternative = getAlternative
    get_probability = getProbability
    get_probabilities = getProbabilities
    is_prior = isPrior
    is_prior_distribution = isPriorDistribution


class MultivariateNormal(ContinuousDistribution):
    """
    Multivariate Normal (Gaussian) Distribution.

    Represents a d-dimensional normal distribution with mean vector mu and
    covariance matrix Sigma. The distribution is characterized by its
    probability density function.

    Can be used standalone or within a Prior for mixture models.

    Args:
        mu (list or array-like): Mean vector (d-dimensional)
        sigma (list or 2D array-like): Covariance matrix (d×d, positive definite)

    Attributes:
        dimension (int): Dimensionality of the distribution

    Examples:
        >>> # 2D normal distribution
        >>> mvn = MultivariateNormal([1.0, 2.0], [[1.0, 0.5], [0.5, 1.0]])
        >>> mvn.getDimension()
        2
        >>> samples = mvn.sample(100)  # Returns 100×2 numpy array
        >>> pdf_value = mvn.evalPDF([1.0, 2.0])
        >>>
        >>> # Use in Prior mixture
        >>> mvn1 = MultivariateNormal([0, 0], [[1, 0], [0, 1]])
        >>> mvn2 = MultivariateNormal([2, 2], [[0.5, 0.2], [0.2, 0.5]])
        >>> prior = Prior([mvn1, mvn2], [0.6, 0.4])
    """

    def __init__(self, *args):
        super().__init__()

        if len(args) == 1:
            # Existing Java object
            self.obj = args[0]
        elif len(args) == 2:
            mu = args[0]
            sigma = args[1]

            # Convert Python lists/arrays to Java Matrix objects
            # CRITICAL: Hide Matrix class from user - use jlineMatrixFromArray
            mu_matrix = jlineMatrixFromArray(mu)
            sigma_matrix = jlineMatrixFromArray(sigma)

            # Create Java MultivariateNormal object
            self.obj = jpype.JPackage('jline').lang.processes.MultivariateNormal(
                mu_matrix, sigma_matrix
            )
        else:
            raise ValueError("MultivariateNormal requires 1 or 2 arguments")

    def getDimension(self):
        """
        Get the dimensionality of the distribution.

        Returns:
            int: Number of dimensions
        """
        return int(self.obj.getDimension())

    def getMeanVector(self):
        """
        Get the mean vector.

        Returns:
            numpy.ndarray: Mean vector (d-dimensional)
        """
        return jlineMatrixToArray(self.obj.getMeanVector()).flatten()

    def getCovariance(self):
        """
        Get the covariance matrix.

        Returns:
            numpy.ndarray: Covariance matrix (d×d)
        """
        return jlineMatrixToArray(self.obj.getCovariance())

    def getCorrelation(self):
        """
        Get the correlation matrix.

        Returns:
            numpy.ndarray: Correlation matrix (d×d)
        """
        return jlineMatrixToArray(self.obj.getCorrelation())

    def sample(self, n=1):
        """
        Generate random samples from the distribution.

        Args:
            n (int): Number of samples to generate (default: 1)

        Returns:
            numpy.ndarray: n×d array of samples
        """
        RandomManager = jpype.JPackage('jline').util.RandomManager
        samples_matrix = self.obj.sampleMatrix(n, RandomManager.getThreadRandomAsRandom())
        return jlineMatrixToArray(samples_matrix)

    def evalPDF(self, x):
        """
        Evaluate probability density function at point(s) x.

        Args:
            x (array-like): Point(s) at which to evaluate PDF
                           Can be 1D (single point) or 2D (multiple points)

        Returns:
            float or numpy.ndarray: PDF value(s)
        """
        x_array = np.array(x)

        if x_array.ndim == 1:
            # Single point - convert to column vector
            x_matrix = jlineMatrixFromArray(x_array.reshape(-1, 1))
            return float(self.obj.evalPDF(x_matrix))
        else:
            # Multiple points - evaluate each row
            pdf_values = []
            for i in range(x_array.shape[0]):
                x_matrix = jlineMatrixFromArray(x_array[i].reshape(-1, 1))
                pdf_values.append(self.obj.evalPDF(x_matrix))
            return np.array(pdf_values)

    def getMarginal(self, indices):
        """
        Extract marginal distribution for subset of dimensions.

        Args:
            indices (list or array): List of dimension indices to keep (0-based)

        Returns:
            MultivariateNormal: Marginal distribution
        """
        indices_java = jpype.JArray(jpype.JInt)(indices)
        return MultivariateNormal(self.obj.getMarginal(indices_java))

    def getMarginalUniv(self, index):
        """
        Extract univariate marginal distribution.

        Args:
            index (int): Dimension index (0-based)

        Returns:
            Normal: Univariate normal distribution for that dimension
        """
        return Normal(self.obj.getMarginalUniv(int(index)))

    # Snake_case aliases for Pythonic interface
    get_dimension = getDimension
    get_mean_vector = getMeanVector
    get_covariance = getCovariance
    get_correlation = getCorrelation
    eval_pdf = evalPDF
    get_marginal = getMarginal
    get_marginal_univ = getMarginalUniv

    @staticmethod
    def fitMeanAndCovariance(mu, sigma):
        """
        Create multivariate normal from mean and covariance.

        Args:
            mu (array-like): Mean vector
            sigma (array-like): Covariance matrix

        Returns:
            MultivariateNormal: Distribution instance
        """
        return MultivariateNormal(mu, sigma)

    fit_mean_and_covariance = fitMeanAndCovariance

