%{
%{
 % @file dist_scale_rate.m
 % @brief Rate-scaled copy of a distribution, preserving its shape.
%}
%}

function scaled = dist_scale_rate(distrib, factor)
%{
%{
 % @brief Returns a new distribution whose rate is FACTOR times the rate of
 %        DISTRIB, i.e. the time-scaled variable X/FACTOR. The scaling is
 %        exact: every moment of order n is divided by FACTOR^n, so the mean
 %        is divided by FACTOR while the SCV, the skewness and the whole
 %        shape of the distribution are preserved.
 %
 %        The scaled object is rebuilt from the parameters of the original,
 %        rather than by rescaling its Markovian representation, so that the
 %        parameter list stays coherent with the distribution family. Solvers
 %        that serialize the model (JMT, LDES) read the parameters, and would
 %        otherwise export the unscaled process.
 %
 %        This is the perturbation primitive of the finite-difference branch
 %        of getSensitivityTable: scaling the rate at a station-class by
 %        (1+h) is exactly the perturbation the derivative d(.)/d(rate) is
 %        taken along.
 %
 % @fn dist_scale_rate(distrib, factor)
 % @param distrib A Distribution object.
 % @param factor  Positive scaling factor for the rate.
 % @return scaled A new Distribution of the same family with rate*FACTOR.
%}
%}

if ~isa(distrib, 'Distribution')
    line_error(mfilename, 'dist_scale_rate expects a Distribution object.');
end
if ~isnumeric(factor) || ~isscalar(factor) || ~isfinite(factor) || factor <= 0
    line_error(mfilename, 'The scaling factor must be a positive finite scalar.');
end

switch class(distrib)
    case 'Exp'
        scaled = Exp(distrib.getParam(1).paramValue * factor);
    case 'Erlang'
        scaled = Erlang(distrib.getParam(1).paramValue * factor, ...
            distrib.getParam(2).paramValue);
    case 'HyperExp'
        p = distrib.getParam(1).paramValue;
        lambda1 = distrib.getParam(2).paramValue;
        lambda2 = distrib.getParam(3).paramValue;
        if isempty(lambda2)
            % n-phase form: param2 holds the whole rate vector.
            scaled = HyperExp(p, lambda1 * factor);
        else
            scaled = HyperExp(p, lambda1 * factor, lambda2 * factor);
        end
    case {'Coxian', 'Cox2'}
        % Completion probabilities are dimensionless and are left untouched;
        % only the phase rates carry the time scale. A Cox2 is a two-phase
        % Coxian and is rebuilt as one, its own constructor being a fit.
        if length(distrib.params) == 3
            scaled = Coxian(distrib.getParam(1).paramValue * factor, ...
                distrib.getParam(2).paramValue * factor, ...
                distrib.getParam(3).paramValue);
        else
            scaled = Coxian(distrib.getParam(1).paramValue * factor, ...
                distrib.getParam(2).paramValue);
        end
    case 'APH'
        scaled = APH(distrib.getParam(1).paramValue, ...
            distrib.getParam(2).paramValue * factor);
    case 'PH'
        scaled = PH(distrib.getParam(1).paramValue, ...
            distrib.getParam(2).paramValue * factor);
    case 'MAP'
        scaled = MAP(distrib.getParam(1).paramValue * factor, ...
            distrib.getParam(2).paramValue * factor);
    case 'MMPP2'
        % Every rate of the modulating chain and of the arrival process is
        % scaled, which time-scales the whole process.
        scaled = MMPP2(distrib.getParam(1).paramValue * factor, ...
            distrib.getParam(2).paramValue * factor, ...
            distrib.getParam(3).paramValue * factor, ...
            distrib.getParam(4).paramValue * factor);
    case 'Det'
        scaled = Det(distrib.getParam(1).paramValue / factor);
    case 'Uniform'
        scaled = Uniform(distrib.getParam(1).paramValue / factor, ...
            distrib.getParam(2).paramValue / factor);
    case 'Gamma'
        % Gamma(shape, scale): the shape is dimensionless.
        scaled = Gamma(distrib.getParam(1).paramValue, ...
            distrib.getParam(2).paramValue / factor);
    case 'Pareto'
        % Pareto(shape, scale): the scale is the minimum of the support.
        scaled = Pareto(distrib.getParam(1).paramValue, ...
            distrib.getParam(2).paramValue / factor);
    case 'Weibull'
        % Weibull params are (1) the scale alpha and (2) the shape r, while
        % the constructor takes (shape, scale).
        scaled = Weibull(distrib.getParam(2).paramValue, ...
            distrib.getParam(1).paramValue / factor);
    case 'Lognormal'
        % X/factor is lognormal with mu - log(factor) and the same sigma.
        scaled = Lognormal(distrib.getParam(1).paramValue - log(factor), ...
            distrib.getParam(2).paramValue);
    case 'NHPP'
        % A piecewise-constant intensity is time-scaled by lambda(t) ->
        % factor*lambda(factor*t), which is the schedule with its rates scaled
        % up and its breakpoints compressed. The time-average rate, which is
        % what sn.rates carries for an NHPP, is then scaled by factor.
        scaled = NHPP(distrib.getBreakpoints() / factor, ...
            distrib.getRates() * factor, distrib.isCyclic());
    case 'Replayer'
        % A trace is scaled sample by sample. A file-backed Replayer is left
        % alone: rescaling it means writing a new trace file, which a
        % sensitivity sweep must not do behind the caller's back.
        if isempty(distrib.data)
            line_error(mfilename, ['Rate scaling is not defined for a ', ...
                'file-backed Replayer: the trace is the parameter. Pass the ', ...
                'samples as an array, Replayer(data), to differentiate it.']);
        end
        scaled = Replayer(distrib.data / factor);
    case 'Immediate'
        scaled = Immediate();
    otherwise
        line_error(mfilename, sprintf(['Rate scaling is not defined for a %s ', ...
            'process. Supported: Exp, Erlang, HyperExp, Coxian, Cox2, APH, PH, ', ...
            'MAP, MMPP2, Det, Uniform, Gamma, Pareto, Weibull, Lognormal, NHPP, ', ...
            'Replayer, Immediate.'], class(distrib)));
end
end
