function entries = line_citations(tokens)
% ENTRIES = LINE_CITATIONS(TOKENS)
%
% Bibliographic references for the algorithms named by TOKENS, as a struct
% array with fields:
%   .key    - bibliography key, as used in doc/latex/biblio.bib
%   .ref    - short reference, author-title-venue-year
%   .covers - one line saying which part of the solution process it covers
%
% TOKENS is a cell array of algorithm or feature names (solver methods such as
% 'bs' or 'comom', transformations such as 'mmt', percentile methods such as
% 'forktail'). Unknown tokens are ignored, so a caller may pass whatever it
% knows about a run. The registry is the manual's method-to-citation table
% (doc/latex/manual.tex) and its bibliography; keep the two in step.
%
% Attribution in LINE is pull-based: nothing is printed during a solve, and a
% user asks for the references with solver.citations() when writing them up.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 1 || isempty(tokens)
    tokens = {};
end
if ~iscell(tokens)
    tokens = {tokens};
end

reg = registry();
entries = struct('key', {}, 'ref', {}, 'covers', {});
seen = {};
for i = 1:numel(tokens)
    t = lower(strtrim(char(tokens{i})));
    if isempty(t)
        continue
    end
    if isKey(reg, t)
        e = reg(t);
    else
        % a family-qualified token falls back to the bare method name
        dot = strfind(t, '.');
        if ~isempty(dot) && isKey(reg, t(dot(1)+1:end))
            e = reg(t(dot(1)+1:end));
        else
            continue
        end
    end
    if any(strcmp(seen, e.key))
        continue    % same paper reached through two tokens
    end
    seen{end+1} = e.key; %#ok<AGROW>
    entries(end+1) = e; %#ok<AGROW>
end
end

function reg = registry()
% Token -> reference. Sources: doc/latex/manual.tex (the solver/method summary
% table and the per-solver method lists) and doc/latex/biblio.bib. A token may
% be qualified by solver family ('nc.mva' is Reiser's convolution, 'mva.exact'
% is Reiser-Lavenberg MVA); the lookup falls back to the bare method name when
% the qualified one is absent.
reg = containers.Map('KeyType','char','ValueType','any');

    function add(tok, key, ref, covers)
        reg(tok) = struct('key', key, 'ref', ref, 'covers', covers);
    end

add('mva.exact', 'reis.lave80', 'M. Reiser, S. Lavenberg, "Mean-Value Analysis of Closed Multichain Queuing Networks", J. ACM 27(2), 1980', 'exact mean queue lengths of the product-form network');
add('exact', 'reis.lave80', 'M. Reiser, S. Lavenberg, "Mean-Value Analysis of Closed Multichain Queuing Networks", J. ACM 27(2), 1980', 'exact mean queue lengths of the product-form network');
add('mva.mva', 'reis.lave80', 'M. Reiser, S. Lavenberg, "Mean-Value Analysis of Closed Multichain Queuing Networks", J. ACM 27(2), 1980', 'exact mean queue lengths of the product-form network');
add('reiser', 'reis.lave80', 'M. Reiser, S. Lavenberg, "Mean-Value Analysis of Closed Multichain Queuing Networks", J. ACM 27(2), 1980', 'exact mean queue lengths of the product-form network');
add('bs', 'Sch79', 'P. J. Schweitzer, "Approximate Analysis of Multiclass Closed Networks of Queues", Int. Conf. Stoch. Control Optim., 1979', 'Bard-Schweitzer fixed point for the mean queue lengths');
add('amva.bs', 'Sch79', 'P. J. Schweitzer, "Approximate Analysis of Multiclass Closed Networks of Queues", Int. Conf. Stoch. Control Optim., 1979', 'Bard-Schweitzer fixed point for the mean queue lengths');
add('aql', 'ZahES88', 'J. Zahorjan, D. L. Eager, H. M. Sweillam, "Accuracy, Speed, and Convergence of Approximate Mean Value Analysis", Perform. Eval. 8, 1988', 'aggregate queue-length AMVA iteration');
add('amva.aql', 'ZahES88', 'J. Zahorjan, D. L. Eager, H. M. Sweillam, "Accuracy, Speed, and Convergence of Approximate Mean Value Analysis", Perform. Eval. 8, 1988', 'aggregate queue-length AMVA iteration');
add('lin', 'ChaN82', 'K. M. Chandy, D. Neuse, "Linearizer: A Heuristic Algorithm for Queuing Network Models of Computing Systems", Commun. ACM 25(2), 1982', 'Linearizer correction of the arrival-instant queue lengths');
add('amva.lin', 'ChaN82', 'K. M. Chandy, D. Neuse, "Linearizer: A Heuristic Algorithm for Queuing Network Models of Computing Systems", Commun. ACM 25(2), 1982', 'Linearizer correction of the arrival-instant queue lengths');
add('gflin', 'ChaN82', 'K. M. Chandy, D. Neuse, "Linearizer: A Heuristic Algorithm for Queuing Network Models of Computing Systems", Commun. ACM 25(2), 1982', 'Linearizer correction of the arrival-instant queue lengths');
add('egflin', 'ChaN82', 'K. M. Chandy, D. Neuse, "Linearizer: A Heuristic Algorithm for Queuing Network Models of Computing Systems", Commun. ACM 25(2), 1982', 'Linearizer correction of the arrival-instant queue lengths');
add('dmlin', 'SilvaM90', 'E. de Souza e Silva, R. R. Muntz, "A Note on the Computational Cost of the Linearizer Algorithm for Queueing Networks", IEEE TC 39(6), 1990', 'de Souza e Silva-Muntz cost reduction of Linearizer');
add('qd', 'casale2015qdamva', 'G. Casale, J. F. Perez, W. Wang, "QD-AMVA: Evaluating Systems with Queue-Dependent Service Requirements", IFIP PERFORMANCE, 2015', 'queue-dependent AMVA for load-dependent stations');
add('amva.qd', 'casale2015qdamva', 'G. Casale, J. F. Perez, W. Wang, "QD-AMVA: Evaluating Systems with Queue-Dependent Service Requirements", IFIP PERFORMANCE, 2015', 'queue-dependent AMVA for load-dependent stations');
add('qdlin', 'casale2015qdamva', 'G. Casale, J. F. Perez, W. Wang, "QD-AMVA: Evaluating Systems with Queue-Dependent Service Requirements", IFIP PERFORMANCE, 2015', 'queue-dependent AMVA for load-dependent stations');
add('softmin', 'casale2015qdamva', 'G. Casale, J. F. Perez, W. Wang, "QD-AMVA: Evaluating Systems with Queue-Dependent Service Requirements", IFIP PERFORMANCE, 2015', 'queue-dependent AMVA for load-dependent stations');
add('qd.oi', 'casale2015qdamva', 'G. Casale, J. F. Perez, W. Wang, "QD-AMVA: Evaluating Systems with Queue-Dependent Service Requirements", IFIP PERFORMANCE, 2015', 'first-order Taylor and Schweitzer steps of the queue-dependent AMVA specialized to order-independent stations');
add('amva.qdoi', 'casale2015qdamva', 'G. Casale, J. F. Perez, W. Wang, "QD-AMVA: Evaluating Systems with Queue-Dependent Service Requirements", IFIP PERFORMANCE, 2015', 'first-order Taylor and Schweitzer steps of the queue-dependent AMVA specialized to order-independent stations');
add('oi', 'BonP03', 'T. Bonald, A. Proutiere, "Insensitive bandwidth sharing in data networks", Queueing Systems 44(1), 2003', 'balanced-fairness balance function of the order-independent station');
add('balancedfairness', 'BonP03', 'T. Bonald, A. Proutiere, "Insensitive bandwidth sharing in data networks", Queueing Systems 44(1), 2003', 'balanced-fairness balance function of the order-independent station');
add('aumannshapley', 'BilH82', 'L. J. Billera, D. C. Heath, "Allocation of Shared Costs: A Set of Axioms Yielding a Unique Procedure", Mathematics of Operations Research 7(1), 1982', 'ray decomposition splitting the order-independent station rate among the job classes');
add('fli', 'WangS00', 'H. Wang, K. C. Sevcik, "Experiments with improved approximate mean value analysis algorithms", Perform. Eval. 39, 2000', 'improved AMVA arrival-instant estimators');
add('amva.fli', 'WangS00', 'H. Wang, K. C. Sevcik, "Experiments with improved approximate mean value analysis algorithms", Perform. Eval. 39, 2000', 'improved AMVA arrival-instant estimators');
add('qli', 'WangS00', 'H. Wang, K. C. Sevcik, "Experiments with improved approximate mean value analysis algorithms", Perform. Eval. 39, 2000', 'improved AMVA arrival-instant estimators');
add('amva.qli', 'WangS00', 'H. Wang, K. C. Sevcik, "Experiments with improved approximate mean value analysis algorithms", Perform. Eval. 39, 2000', 'improved AMVA arrival-instant estimators');
add('conway', 'Con89', 'A. E. Conway, "Fast Approximate Solution of Queueing Networks with Multi-Server Chain-Dependent FCFS Queues", 1989', 'multi-server chain-dependent AMVA');
add('linearizerms', 'Con89', 'A. E. Conway, "Fast Approximate Solution of Queueing Networks with Multi-Server Chain-Dependent FCFS Queues", 1989', 'multi-server chain-dependent AMVA');
add('suri', 'suri2007approximate', 'R. Suri, S. K. Sahu, M. Vernon, "Approximate Mean Value Analysis for Closed Queuing Networks with Multiple-Server Stations", IERC, 2007', 'multi-server station correction in AMVA');
add('schmidt', 'suri2007approximate', 'R. Suri, S. K. Sahu, M. Vernon, "Approximate Mean Value Analysis for Closed Queuing Networks with Multiple-Server Stations", IERC, 2007', 'multi-server station correction in AMVA');
add('schmidt-ext', 'suri2007approximate', 'R. Suri, S. K. Sahu, M. Vernon, "Approximate Mean Value Analysis for Closed Queuing Networks with Multiple-Server Stations", IERC, 2007', 'multi-server station correction in AMVA');
add('schmidtext', 'suri2007approximate', 'R. Suri, S. K. Sahu, M. Vernon, "Approximate Mean Value Analysis for Closed Queuing Networks with Multiple-Server Stations", IERC, 2007', 'multi-server station correction in AMVA');
add('seidmann', 'seidmann1987computerized', 'A. Seidmann, P. J. Schweitzer, S. Shalev-Oren, "Computerized closed queueing network models of flexible manufacturing systems", Large Scale Systems 12, 1987', 'multi-server to single-server flow-equivalent reduction');
add('cl', 'EagL88', 'D. L. Eager, J. N. Lipscomb, "The AMVA priority approximation", Perform. Eval. 8, 1988', 'priority approximation inside AMVA');
add('chandy-lakshmi', 'EagL88', 'D. L. Eager, J. N. Lipscomb, "The AMVA priority approximation", Perform. Eval. 8, 1988', 'priority approximation inside AMVA');
add('amva.cl', 'EagL88', 'D. L. Eager, J. N. Lipscomb, "The AMVA priority approximation", Perform. Eval. 8, 1988', 'priority approximation inside AMVA');
add('shadow', 'Sev77', 'K. Sevcik, "Priority Scheduling Disciplines in Queuing Network Models of Computer Systems", IFIP Congress, 1977', 'shadow-server treatment of priority scheduling');
add('zhou', 'Woo22', 'S. Zhou, M. Woodside, "A Multiserver Approximation for Cloud Scaling Analysis", ICPE Companion, 2022', 'multiserver scaling approximation');
add('qna', 'whitt1983qna', 'W. Whitt, "The Queueing Network Analyzer", Bell Syst. Tech. J. 62, 1983', 'two-moment decomposition of the open network');
add('rqna', 'whittyou2018rqna', 'W. Whitt, W. You, "A Robust Queueing Network Analyzer Based on Indices of Dispersion", Naval Research Logistics 69, 2022', 'index-of-dispersion decomposition and robust queueing bounds');
add('highvar', 'BonW86', 'A. B. Bondi, W. Whitt, "The influence of service-time variability in a closed network of queues", Perform. Eval. 6, 1986', 'high-variability service correction of the demands');
add('interp', 'BonW86', 'A. B. Bondi, W. Whitt, "The influence of service-time variability in a closed network of queues", Perform. Eval. 6, 1986', 'high-variability service correction of the demands');
add('kraemer', 'KraLB78', 'W. Kraemer, M. Langenbach-Belz, "Approximate Formulae for General Single Server Systems with Single and Batch Arrivals", Angewandte Informatik 9, 1978', 'G/G/1 waiting time approximation');
add('klb', 'KraLB78', 'W. Kraemer, M. Langenbach-Belz, "Approximate Formulae for General Single Server Systems with Single and Batch Arrivals", Angewandte Informatik 9, 1978', 'G/G/1 waiting time approximation');
add('nc.ca', 'Cas09', 'G. Casale, "CoMoM: Efficient Class-Oriented Evaluation of Multiclass Performance Models", IEEE TSE 35(2), 2009', 'class-oriented recursion for the normalizing constant');
add('ca', 'Cas09', 'G. Casale, "CoMoM: Efficient Class-Oriented Evaluation of Multiclass Performance Models", IEEE TSE 35(2), 2009', 'class-oriented recursion for the normalizing constant');
% Loss networks with finite capacity regions (SolverNC, solver_nc_lossn_analyzer)
add('lossn.exact', 'ManSik07', 'D. Manjunath, B. Sikdar, "Integral Expressions for the Numerical Evaluation of Product Form Expressions Over Irregular Multidimensional Integer Spaces"', 'exact normalizing constant over a state space cut by linear integer constraints');
add('lossn.ms', 'ManSik07', 'D. Manjunath, B. Sikdar, "Integral Expressions for the Numerical Evaluation of Product Form Expressions Over Irregular Multidimensional Integer Spaces"', 'exact normalizing constant over a state space cut by linear integer constraints');
add('lossn.erlangfp', 'Kelly91', 'F. P. Kelly, "Loss Networks", Ann. Appl. Probab. 1(3), 1991', 'reduced-load (Erlang fixed point) approximation of link blocking');
add('erlangfp', 'Kelly91', 'F. P. Kelly, "Loss Networks", Ann. Appl. Probab. 1(3), 1991', 'reduced-load (Erlang fixed point) approximation of link blocking');
add('lossn.mci', 'RosWan92', 'K. W. Ross, J. Wang, "Monte Carlo Summation Applied to Product-Form Loss Networks", Prob. Eng. Inf. Sci. 6, 1992', 'Monte Carlo summation of the loss-network normalizing constant');
add('mci', 'RosWan92', 'K. W. Ross, J. Wang, "Monte Carlo Summation Applied to Product-Form Loss Networks", Prob. Eng. Inf. Sci. 6, 1992', 'Monte Carlo summation of the loss-network normalizing constant');

add('comom', 'Cas09', 'G. Casale, "CoMoM: Efficient Class-Oriented Evaluation of Multiclass Performance Models", IEEE TSE 35(2), 2009', 'class-oriented recursion for the normalizing constant');
add('nc.comom', 'Cas09', 'G. Casale, "CoMoM: Efficient Class-Oriented Evaluation of Multiclass Performance Models", IEEE TSE 35(2), 2009', 'class-oriented recursion for the normalizing constant');
add('clw', 'ChoLW95', 'G. L. Choudhury, K. K. Leung, W. Whitt, "Calculating Normalization Constants of Closed Queuing Networks by Numerically Inverting Their Generating Functions", J. ACM 42, 1995', 'numerical inversion of the generating function');
add('nc.clw', 'ChoLW95', 'G. L. Choudhury, K. K. Leung, W. Whitt, "Calculating Normalization Constants of Closed Queuing Networks by Numerically Inverting Their Generating Functions", J. ACM 42, 1995', 'numerical inversion of the generating function');
add('nc.mva', 'Rei81', 'M. Reiser, "Mean-Value Analysis and Convolution Method for Queue-Dependent Servers in Closed Queueing Networks", Perform. Eval. 1, 1981', 'convolution and MVA with queue-dependent servers');
add('le', 'Cas17', 'G. Casale, "Accelerating Performance Inference over Closed Systems by Asymptotic Methods", ACM SIGMETRICS, 2017', 'asymptotic expansion of the normalizing constant integral');
add('ls', 'Cas17', 'G. Casale, "Accelerating Performance Inference over Closed Systems by Asymptotic Methods", ACM SIGMETRICS, 2017', 'asymptotic expansion of the normalizing constant integral');
add('cub', 'Cas17', 'G. Casale, "Accelerating Performance Inference over Closed Systems by Asymptotic Methods", ACM SIGMETRICS, 2017', 'asymptotic expansion of the normalizing constant integral');
add('nc.le', 'Cas17', 'G. Casale, "Accelerating Performance Inference over Closed Systems by Asymptotic Methods", ACM SIGMETRICS, 2017', 'asymptotic expansion of the normalizing constant integral');
add('nc.ls', 'Cas17', 'G. Casale, "Accelerating Performance Inference over Closed Systems by Asymptotic Methods", ACM SIGMETRICS, 2017', 'asymptotic expansion of the normalizing constant integral');
add('nc.cub', 'Cas17', 'G. Casale, "Accelerating Performance Inference over Closed Systems by Asymptotic Methods", ACM SIGMETRICS, 2017', 'asymptotic expansion of the normalizing constant integral');
add('rd', 'CasHH21', 'G. Casale, P. G. Harrison, O. W. Hong, "Facilitating Load-Dependent Queueing Analysis Through Factorization", Perform. Eval., 2021', 'factorization of the load-dependent normalizing constant');
add('nrp', 'CasHH21', 'G. Casale, P. G. Harrison, O. W. Hong, "Facilitating Load-Dependent Queueing Analysis Through Factorization", Perform. Eval., 2021', 'factorization of the load-dependent normalizing constant');
add('nrl', 'CasHH21', 'G. Casale, P. G. Harrison, O. W. Hong, "Facilitating Load-Dependent Queueing Analysis Through Factorization", Perform. Eval., 2021', 'factorization of the load-dependent normalizing constant');
add('comomld', 'CasHH21', 'G. Casale, P. G. Harrison, O. W. Hong, "Facilitating Load-Dependent Queueing Analysis Through Factorization", Perform. Eval., 2021', 'factorization of the load-dependent normalizing constant');
add('kt', 'KneT92', 'C. Knessl, C. Tier, "Asymptotic Expansions for Large Closed Queueing Networks with Multiple Job Classes", IEEE TC 41(4), 1992', 'asymptotic expansion for large populations');
add('panacea', 'McKM84', 'J. McKenna, D. Mitra, "Asymptotic Expansions and Integral Representations of Moments of Queue Lengths in Closed Markovian Networks", J. ACM 31, 1984', 'integral representation of the queue-length moments');
add('panaceald', 'MitM86', 'D. Mitra, J. McKenna, "Asymptotic Expansions for Closed Markovian Networks with State-Dependent Service Rates", J. ACM 33(3), 1986', 'load-dependent PANACEA expansion and its pseudonetwork coefficients');
add('psrespt', 'MitMo83', 'D. Mitra, J. A. Morrison, "Asymptotic Expansions of Moments of the Waiting Time in Closed and Open Processor-Sharing Systems with Multiple Job Classes", Adv. Appl. Prob. 15(4), 1983', 'sojourn-time moments at a multiclass processor-sharing station');
add('mm1ps', 'MitMo83', 'D. Mitra, J. A. Morrison, "Asymptotic Expansions of Moments of the Waiting Time in Closed and Open Processor-Sharing Systems with Multiple Job Classes", Adv. Appl. Prob. 15(4), 1983', 'sojourn-time moments at a multiclass processor-sharing station');
add('mem', 'Kou94', 'D. D. Kouvatsos, "Entropy Maximisation and Queueing Network Models", Annals of Operations Research 48, 1994', 'maximum-entropy approximation of the network');
add('mem.blocking', 'TahMB99', 'H. Tahilramani, D. Manjunath, S. K. Bose, "Approximate Analysis of Open Network of GE/GE/m/N Queues with Transfer Blocking", MASCOTS, 1999', 'holding-node expansion that makes transfer blocking work conserving');
add('recal', 'ConG86', 'A. E. Conway, N. D. Georganas, "RECAL: A New Efficient Algorithm for the Exact Analysis of Multiple-Chain Closed Queueing Networks", J. ACM 33, 1986', 'recursive exact evaluation by chain');
add('mvac', 'CSL89', 'A. E. Conway, E. de Souza e Silva, S. S. Lavenberg, "Mean Value Analysis by Chain of Product Form Queueing Networks", IEEE Trans. Computers 38(3), 1989', 'exact mean value analysis by chain (SolverMVA method mvac)');
add('conv', 'Sau83', 'C. H. Sauer, "Computational Algorithms for State-Dependent Queueing Networks", ACM TOCS 1(1), 1983', 'convolution with chain-dependent service rates');
add('nc.conv', 'Sau83', 'C. H. Sauer, "Computational Algorithms for State-Dependent Queueing Networks", ACM TOCS 1(1), 1983', 'convolution with chain-dependent service rates');
add('aba', 'BolGMT06', 'G. Bolch, S. Greiner, H. de Meer, K. S. Trivedi, "Queueing Networks and Markov Chains", Wiley, 2006', 'asymptotic bounds on throughput and response time');
add('aba.upper', 'BolGMT06', 'G. Bolch, S. Greiner, H. de Meer, K. S. Trivedi, "Queueing Networks and Markov Chains", Wiley, 2006', 'asymptotic bounds on throughput and response time');
add('aba.lower', 'BolGMT06', 'G. Bolch, S. Greiner, H. de Meer, K. S. Trivedi, "Queueing Networks and Markov Chains", Wiley, 2006', 'asymptotic bounds on throughput and response time');
add('bjb', 'CasMS08', 'G. Casale, R. R. Muntz, G. Serazzi, "Geometric Bounds: A Noniterative Analysis Technique for Closed Queueing Networks", IEEE TC 57(6), 2008', 'geometric and balanced-job bounds');
add('gb', 'CasMS08', 'G. Casale, R. R. Muntz, G. Serazzi, "Geometric Bounds: A Noniterative Analysis Technique for Closed Queueing Networks", IEEE TC 57(6), 2008', 'geometric and balanced-job bounds');
add('pb', 'CasMS08', 'G. Casale, R. R. Muntz, G. Serazzi, "Geometric Bounds: A Noniterative Analysis Technique for Closed Queueing Networks", IEEE TC 57(6), 2008', 'geometric and balanced-job bounds');
add('bjb.upper', 'CasMS08', 'G. Casale, R. R. Muntz, G. Serazzi, "Geometric Bounds: A Noniterative Analysis Technique for Closed Queueing Networks", IEEE TC 57(6), 2008', 'geometric and balanced-job bounds');
add('bjb.lower', 'CasMS08', 'G. Casale, R. R. Muntz, G. Serazzi, "Geometric Bounds: A Noniterative Analysis Technique for Closed Queueing Networks", IEEE TC 57(6), 2008', 'geometric and balanced-job bounds');
add('gb.upper', 'CasMS08', 'G. Casale, R. R. Muntz, G. Serazzi, "Geometric Bounds: A Noniterative Analysis Technique for Closed Queueing Networks", IEEE TC 57(6), 2008', 'geometric and balanced-job bounds');
add('gb.lower', 'CasMS08', 'G. Casale, R. R. Muntz, G. Serazzi, "Geometric Bounds: A Noniterative Analysis Technique for Closed Queueing Networks", IEEE TC 57(6), 2008', 'geometric and balanced-job bounds');
add('pb.upper', 'CasMS08', 'G. Casale, R. R. Muntz, G. Serazzi, "Geometric Bounds: A Noniterative Analysis Technique for Closed Queueing Networks", IEEE TC 57(6), 2008', 'geometric and balanced-job bounds');
add('pb.lower', 'CasMS08', 'G. Casale, R. R. Muntz, G. Serazzi, "Geometric Bounds: A Noniterative Analysis Technique for Closed Queueing Networks", IEEE TC 57(6), 2008', 'geometric and balanced-job bounds');
add('sb', 'Harel1999', 'A. Harel, S. Namn, J. Sturm, "Simple bounds for closed queueing networks", Queueing Systems 31, 1999', 'simple closed-network bounds');
add('sb.upper', 'Harel1999', 'A. Harel, S. Namn, J. Sturm, "Simple bounds for closed queueing networks", Queueing Systems 31, 1999', 'simple closed-network bounds');
add('sb.lower', 'Harel1999', 'A. Harel, S. Namn, J. Sturm, "Simple bounds for closed queueing networks", Queueing Systems 31, 1999', 'simple closed-network bounds');
add('mwba', 'MajW98', 'S. Majumdar, C. M. Woodside, "Robust bounds and throughput guarantees for closed multiclass queueing networks", Perform. Eval. 32, 1998', 'robust throughput bounds');
add('mwba.upper', 'MajW98', 'S. Majumdar, C. M. Woodside, "Robust bounds and throughput guarantees for closed multiclass queueing networks", Perform. Eval. 32, 1998', 'robust throughput bounds');
add('mwba.lower', 'MajW98', 'S. Majumdar, C. M. Woodside, "Robust bounds and throughput guarantees for closed multiclass queueing networks", Perform. Eval. 32, 1998', 'robust throughput bounds');
add('balanced', 'LazZGS84', 'E. D. Lazowska, J. Zahorjan, G. S. Graham, K. C. Sevcik, "Quantitative System Performance", Prentice-Hall, 1984', 'operational analysis and balanced-system bounds');
add('mna', 'ZhuC24', 'Z. Li, G. Casale, "Matrix Network Analyzer: A New Decomposition Algorithm for Phase-type Queueing Networks", ICPE Companion, 2024', 'phase-type network decomposition');
add('inap', 'CasH13', 'G. Casale, P. G. Harrison, "AutoCAT: Automated Product-Form Solution of Stochastic Models", MAM in Stochastic Models 27, 2013', 'RCAT product-form solution of the cooperating processes');
add('inapplus', 'CasH13', 'G. Casale, P. G. Harrison, "AutoCAT: Automated Product-Form Solution of Stochastic Models", MAM in Stochastic Models 27, 2013', 'RCAT product-form solution of the cooperating processes');
add('inapinf', 'MarinRB12', 'A. Marin, S. Rota Bulo, S. Balsamo, "A numerical algorithm for the decomposition of cooperating structured Markov processes", IEEE MASCOTS, 2012', 'matrix-geometric decomposition on the infinite state space');
add('qbd', 'Hor17', 'G. Horvath, M. Telek, "BuTools 2: A Rich Toolbox for Markovian Performance Evaluation", VALUETOOLS, 2017', 'quasi-birth-death and matrix-analytic routines');
add('mam.qbd', 'Hor17', 'G. Horvath, M. Telek, "BuTools 2: A Rich Toolbox for Markovian Performance Evaluation", VALUETOOLS, 2017', 'quasi-birth-death and matrix-analytic routines');
add('mg1.fb', 'WieH03', 'A. Wierman, M. Harchol-Balter, "Classifying scheduling policies with respect to unfairness in an M/GI/1", ACM SIGMETRICS, 2003', 'size-based M/G/1 scheduling response times');
add('mg1.lrpt', 'WieH03', 'A. Wierman, M. Harchol-Balter, "Classifying scheduling policies with respect to unfairness in an M/GI/1", ACM SIGMETRICS, 2003', 'size-based M/G/1 scheduling response times');
add('mg1.psjf', 'WieH03', 'A. Wierman, M. Harchol-Balter, "Classifying scheduling policies with respect to unfairness in an M/GI/1", ACM SIGMETRICS, 2003', 'size-based M/G/1 scheduling response times');
add('mg1.srpt', 'WieH03', 'A. Wierman, M. Harchol-Balter, "Classifying scheduling policies with respect to unfairness in an M/GI/1", ACM SIGMETRICS, 2003', 'size-based M/G/1 scheduling response times');
add('fb', 'WieH03', 'A. Wierman, M. Harchol-Balter, "Classifying scheduling policies with respect to unfairness in an M/GI/1", ACM SIGMETRICS, 2003', 'size-based M/G/1 scheduling response times');
add('lrpt', 'WieH03', 'A. Wierman, M. Harchol-Balter, "Classifying scheduling policies with respect to unfairness in an M/GI/1", ACM SIGMETRICS, 2003', 'size-based M/G/1 scheduling response times');
add('psjf', 'WieH03', 'A. Wierman, M. Harchol-Balter, "Classifying scheduling policies with respect to unfairness in an M/GI/1", ACM SIGMETRICS, 2003', 'size-based M/G/1 scheduling response times');
add('srpt', 'WieH03', 'A. Wierman, M. Harchol-Balter, "Classifying scheduling policies with respect to unfairness in an M/GI/1", ACM SIGMETRICS, 2003', 'size-based M/G/1 scheduling response times');
add('mg1.setf', 'NuyW08', 'M. Nuyens, A. Wierman, "The Foreground-Background queue: A survey", Perform. Eval. 65, 2008', 'foreground-background M/G/1 response times');
add('setf', 'NuyW08', 'M. Nuyens, A. Wierman, "The Foreground-Background queue: A survey", Perform. Eval. 65, 2008', 'foreground-background M/G/1 response times');
add('mapm1ps', 'MasuyamaTakine2003', 'H. Masuyama, T. Takine, "Sojourn time distribution in a MAP/M/1 processor-sharing queue", Oper. Res. Lett. 31, 2003', 'MAP/M/1-PS sojourn time distribution');
add('mmt', 'dobre24', 'R.-A. Dobre, Z. Niu, G. Casale, "Approximating Fork-Join Systems via Mixed Model Transformations", ICPE Companion, 2024', 'mixed model transformation of the fork-join network');
add('fjt', 'dobre24', 'R.-A. Dobre, Z. Niu, G. Casale, "Approximating Fork-Join Systems via Mixed Model Transformations", ICPE Companion, 2024', 'mixed model transformation of the fork-join network');
add('ht', 'heidelberger1982queueing', 'P. Heidelberger, K. Trivedi, "Queueing network models for parallel processing with asynchronous tasks", IEEE TC C-31(11), 1982', 'fork-join transformation with auxiliary asynchronous tasks');
add('heidelberger-trivedi', 'heidelberger1982queueing', 'P. Heidelberger, K. Trivedi, "Queueing network models for parallel processing with asynchronous tasks", IEEE TC C-31(11), 1982', 'fork-join transformation with auxiliary asynchronous tasks');
add('forktail', 'NguALCJ18', 'M. Nguyen, S. Alesawi, N. Li, H. Che, H. Jiang, "ForkTail: A Black-Box Fork-Join Tail Latency Prediction Model", ACM HPDC, 2018', 'response time tail of the fork-join request from the branch moments');
add('qiu', 'QiuPH15', 'Z. Qiu, J. F. Perez, P. G. Harrison, "Beyond the Mean in Fork-Join Queues: Efficient Approximation for Response-Time Tails", IFIP PERFORMANCE, 2015', 'response time tail of a homogeneous fork-join network');
add('ctmc', 'BolGMT06', 'G. Bolch, S. Greiner, H. de Meer, K. S. Trivedi, "Queueing Networks and Markov Chains", Wiley, 2006', 'CTMC formulation, uniformization and its stationary solution');
add('uniformization', 'BolGMT06', 'G. Bolch, S. Greiner, H. de Meer, K. S. Trivedi, "Queueing Networks and Markov Chains", Wiley, 2006', 'CTMC formulation, uniformization and its stationary solution');
add('courtois', 'courtois1977decomposability', 'P. J. Courtois, "Decomposability: Queueing and Computer System Applications", Academic Press, 1977', 'nearly-completely-decomposable aggregation of the chain');
add('kms', 'koury1984iterative', 'J. R. Koury, D. F. McAllister, W. J. Stewart, "Iterative methods for computing stationary distributions of nearly completely decomposable Markov chains", SIAM J. Alg. Disc. Meth. 5, 1984', 'iterative aggregation-disaggregation of the chain');
add('takahashi', 'takahashi1975iterative', 'Y. Takahashi, "A lumping method for numerical calculations of stationary distributions of Markov chains", 1975', 'lumping-based iterative solution of the chain');
add('qrf', 'CasNPS16', 'G. Casale, V. De Nitto Persone, E. Smirni, "QRF: An Optimization-Based Framework for Evaluating Complex Stochastic Networks", ACM TOMACS 26, 2016', 'quadratic reduction bounds on the chain');
add('qrf.mmi', 'CasNPS16', 'G. Casale, V. De Nitto Persone, E. Smirni, "QRF: An Optimization-Based Framework for Evaluating Complex Stochastic Networks", ACM TOMACS 26, 2016', 'quadratic reduction bounds on the chain');
add('qrf.mem', 'CasNPS16', 'G. Casale, V. De Nitto Persone, E. Smirni, "QRF: An Optimization-Based Framework for Evaluating Complex Stochastic Networks", ACM TOMACS 26, 2016', 'quadratic reduction bounds on the chain');
add('qrf.bas', 'CasNPS16', 'G. Casale, V. De Nitto Persone, E. Smirni, "QRF: An Optimization-Based Framework for Evaluating Complex Stochastic Networks", ACM TOMACS 26, 2016', 'quadratic reduction bounds on the chain');
add('qrf.rsrd', 'CasNPS16', 'G. Casale, V. De Nitto Persone, E. Smirni, "QRF: An Optimization-Based Framework for Evaluating Complex Stochastic Networks", ACM TOMACS 26, 2016', 'quadratic reduction bounds on the chain');
add('ssa', 'Gill77', 'D. T. Gillespie, "Exact stochastic simulation of coupled chemical reactions", J. Phys. Chem. 81(25), 1977', 'stochastic simulation of the Markov process sample path');
add('firingdep', 'mars.ea84', 'M. Ajmone Marsan, G. Conte, G. Balbo, "A class of generalized stochastic Petri nets for the performance evaluation of multiprocessor systems", ACM TOCS 2(2), 1984', 'marking-dependent transition firing rates (setFiringRateDependence)');
add('oi', 'BonP03', 'T. Bonald, A. Proutiere, "Insensitive bandwidth sharing in data networks", Queueing Systems 44(1), 2003', 'balanced-fairness balance function of the order-independent station');
add('nc.oi', 'BonP03', 'T. Bonald, A. Proutiere, "Insensitive bandwidth sharing in data networks", Queueing Systems 44(1), 2003', 'balanced-fairness balance function of the order-independent station');
add('mva.oi', 'BonP03', 'T. Bonald, A. Proutiere, "Insensitive bandwidth sharing in data networks", Queueing Systems 44(1), 2003', 'balanced-fairness balance function of the order-independent station');
add('pas', 'ComD21', 'C. Comte, J.-P. Dorsman, "Pass-and-swap queues", Queueing Systems, 2021 (arXiv:2009.12299)', 'per-communicating-class product form of the pass-and-swap network');
add('nrm', 'And07', 'D. F. Anderson, "A modified next reaction method for simulating chemical systems with time dependent propensities and delays", J. Chem. Phys. 127, 2007', 'next-reaction method for the simulation');
add('fld', 'PerC17', 'J. F. Perez, G. Casale, "LINE: Evaluating Software Applications in Unreliable Environments", IEEE Trans. Reliability 66(3), 2017', 'mean-field fluid ODEs for the queueing network');
add('fluid', 'PerC17', 'J. F. Perez, G. Casale, "LINE: Evaluating Software Applications in Unreliable Environments", IEEE Trans. Reliability 66(3), 2017', 'mean-field fluid ODEs for the queueing network');
add('statedep', 'PerC17', 'J. F. Perez, G. Casale, "LINE: Evaluating Software Applications in Unreliable Environments", IEEE Trans. Reliability 66(3), 2017', 'mean-field fluid ODEs for the queueing network');
add('closing', 'PerC17', 'J. F. Perez, G. Casale, "LINE: Evaluating Software Applications in Unreliable Environments", IEEE Trans. Reliability 66(3), 2017', 'mean-field fluid ODEs for the queueing network');
add('matrix', 'RuuskanenBAC21', 'J. Ruuskanen, T. Berner, K.-E. Arzen, A. Cervin, "Improving the mean-field fluid model of processor sharing queueing networks", Perform. Eval. 151, 2021', 'matrix form of the processor-sharing fluid model');
add('rmf', 'GastH16', 'N. Gast, B. Van Houdt, "Transient and steady-state regime of a family of list-based cache replacement algorithms", Queueing Syst. 83, 2016', 'refined mean-field approximation');
add('fluid.rmf', 'GastH16', 'N. Gast, B. Van Houdt, "Transient and steady-state regime of a family of list-based cache replacement algorithms", Queueing Syst. 83, 2016', 'refined mean-field approximation');
add('sfifo.rmf', 'GastH16', 'N. Gast, B. Van Houdt, "Transient and steady-state regime of a family of list-based cache replacement algorithms", Queueing Syst. 83, 2016', 'position-resolved mean field for the strict FIFO(m) replacement variant');
add('fifo.rmf.tran', 'GastH16', 'N. Gast, B. Van Houdt, "Transient and steady-state regime of a family of list-based cache replacement algorithms", Queueing Syst. 83, 2016', 'position-resolved mean-field transient for FIFO(m) (steady state equals RANDOM(m))');
add('diffusion', 'BolGMT06', 'G. Bolch, S. Greiner, H. de Meer, K. S. Trivedi, "Queueing Networks and Markov Chains", Wiley, 2006', 'diffusion approximation of the queue-length process');
add('jmt', 'BerCS07', 'M. Bertoli, G. Casale, G. Serazzi, "The JMT Simulator for Performance Evaluation of Non-Product-Form Queueing Networks", ANSS, 2007', 'discrete-event simulation and JMVA analysis');
add('jsim', 'BerCS07', 'M. Bertoli, G. Casale, G. Serazzi, "The JMT Simulator for Performance Evaluation of Non-Product-Form Queueing Networks", ANSS, 2007', 'discrete-event simulation and JMVA analysis');
add('jmva', 'BerCS07', 'M. Bertoli, G. Casale, G. Serazzi, "The JMT Simulator for Performance Evaluation of Non-Product-Form Queueing Networks", ANSS, 2007', 'discrete-event simulation and JMVA analysis');
add('ln', 'roli.sevc95', 'J. A. Rolia, K. C. Sevcik, "The Method of Layers", IEEE TSE 21(8), 1995', 'layer decomposition of the layered queueing network');
add('layers', 'roli.sevc95', 'J. A. Rolia, K. C. Sevcik, "The Method of Layers", IEEE TSE 21(8), 1995', 'layer decomposition of the layered queueing network');
add('lqns', 'lqns12', 'G. Franks, P. Maly, M. Woodside, D. C. Petriu, A. Hubbard, M. Mroz, "Layered Queueing Network Solver and Simulator User Manual", Carleton University, 2012', 'layered queueing solver and simulator');
add('ln.dec', 'fran.ea09', 'G. Franks, T. Al-Omari, M. Woodside, O. Das, S. Derisavi, "Enhanced Modeling and Solution of Layered Queueing Networks", IEEE TSE 35(2), 2009', 'enhanced decomposition of the layers');
add('enhanced', 'fran.ea09', 'G. Franks, T. Al-Omari, M. Woodside, O. Das, S. Derisavi, "Enhanced Modeling and Solution of Layered Queueing Networks", IEEE TSE 35(2), 2009', 'enhanced decomposition of the layers');
add('ln.fluid', 'trib13', 'M. Tribastone, "A Fluid Model for Layered Queueing Networks", IEEE TSE 39(6), 2013', 'fluid model of the layered network');
add('env', 'casa.trib11', 'G. Casale, M. Tribastone, "Fluid Analysis of Queueing in Two-Stage Random Environments", QEST, 2011', 'queueing in a random environment');
add('env.blend', 'pere.casa13', 'J. F. Perez, G. Casale, "Assessing SLA Compliance from Palladio Component Models", MICAS, 2013', 'environment-stage blending of the metrics');
add('tree', 'rice76', 'J. R. Rice, "The Algorithm Selection Problem", Advances in Computers 15, 1976', 'per-instance selection of the solver from model features');
add('auto.tree', 'rice76', 'J. R. Rice, "The Algorithm Selection Problem", Advances in Computers 15, 1976', 'per-instance selection of the solver from model features');
add('cart', 'brei.ea84', 'L. Breiman, J. Friedman, R. Olshen, C. Stone, "Classification and Regression Trees", Wadsworth, 1984', 'decision tree fitted offline to the selection map');
end
