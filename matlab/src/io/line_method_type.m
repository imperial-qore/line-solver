function label = line_method_type(solvername, method)
% LABEL = LINE_METHOD_TYPE(SOLVERNAME, METHOD)
%
% Classification of a solution method, as printed in the solver banner:
%
%   <accuracy>, <randomness>
%
% with ACCURACY in {exact, approximate, bound} and RANDOMNESS in
% {deterministic, randomized}. SOLVERNAME is the banner solver name ('MVA',
% 'NC', ...) and METHOD the resolved method label, which may carry the
% 'default/' prefix produced by the banner ('default/exact').
%
% Conventions, applied uniformly across the four codebases:
%   - exact:       the algorithm targets the metric with no modeling
%                  approximation; numerical truncation and floating-point
%                  error do not make a method approximate. Integral
%                  representations and transform inversions are exact,
%                  asymptotic expansions are not.
%   - approximate: the algorithm introduces a heuristic, an asymptotic
%                  expansion, a decomposition, or a statistical estimate.
%   - bound:       the algorithm returns a formal one-sided bound on the
%                  metric, not a point estimate: the value is guaranteed to
%                  lie on the stated side of the exact one, and the two sides
%                  of a family bracket it. The side is read off the method
%                  label and printed with it ('gb.upper' -> 'upper bound').
%   - randomized:  the algorithm consumes pseudo-random numbers, so two runs
%                  agree only if the seed does.
%   Perfect sampling (cftp) is classified by the law it samples from, which is
%   the stationary one, hence exact; the ordinary simulators are approximate
%   because a finite horizon leaves warm-up bias on top of the sampling error.
%
% Lookup order: '<solver>.<method>', '<method>', the tail after each dot of
% METHOD (longest suffix first), the head before its first dot, '<solver>',
% then the global default 'approximate, deterministic'. The unknown-method
% default is deliberately the conservative one: claiming exactness that a
% method does not have is the costlier error.
%
% Keep this registry in step with line_citations.m and with its twins
% MethodType.java, solvers/base.py:method_type and cpp method_type.h.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if nargin < 1
    solvername = '';
end
if nargin < 2
    method = '';
end
solvername = lower(strtrim(char(solvername)));
solvername = regexprep(solvername, '^solver', '');
method = lower(strtrim(char(method)));

% the banner label is 'default/<resolved>' once a default has been resolved
slash = strfind(method, '/');
if ~isempty(slash)
    method = method(slash(end)+1:end);
end

reg = registry();
cand = {};
if ~isempty(solvername) && ~isempty(method)
    cand{end+1} = [solvername '.' method];
end
if ~isempty(method)
    cand{end+1} = method;
    dot = strfind(method, '.');
    for d = 1:numel(dot)    % 'a.b.c' -> 'b.c', 'c'
        cand{end+1} = method(dot(d)+1:end); %#ok<AGROW>
    end
    if ~isempty(dot)
        cand{end+1} = method(1:dot(1)-1); %#ok<AGROW>
    end
end
if ~isempty(solvername)
    % '#' keeps the per-solver default out of the method namespace: 'mva' is
    % also a method name, and an unknown SolverMVA method must not inherit it
    cand{end+1} = ['#' solvername];
end

for i = 1:numel(cand)
    if isKey(reg, cand{i})
        e = reg{cand{i}};
        label = sprintf('%s, %s', boundSide(e.accuracy, method), e.randomness);
        return
    end
end
label = 'approximate, deterministic';
end

function accuracy = boundSide(accuracy, method)
% A bound is reported with the side it lies on, which the method label carries
% as its last component ('gb.upper' -> 'upper bound'). A family with no sided
% variant (the qrf reductions) stays the unqualified 'bound'.
if strcmp(accuracy, 'bound')
    if endsWith(method, '.upper')
        accuracy = 'upper bound';
    elseif endsWith(method, '.lower')
        accuracy = 'lower bound';
    end
end
end

function reg = registry()
% Method name -> classification. Method names are solver methods as they appear in the
% banner, optionally qualified by solver family ('jmt.jsim'). A family entry
% ('amva', 'dec', 'jmva') classifies every method under it that has no entry
% of its own, which is what keeps a newly added variant correctly labelled.
reg = configureDictionary('string','cell');

    function add(toks, accuracy, randomness)
        if ~iscell(toks)
            toks = {toks};
        end
        for k = 1:numel(toks)
            reg{toks{k}} = struct('accuracy', accuracy, 'randomness', randomness);
        end
    end

% ---- exact, deterministic ------------------------------------------------
% product-form evaluation: MVA, convolution, RECAL/CoMoM and the
% load-dependent factorizations, plus the exact single-queue closed forms
add({'exact','mva','mvac','recal','conv','ca','comom','comomld', ...
    'rd','nrp','nrl','nre','clw','gleint','mmint2','lcfsqn.ca','rgf','dnc','ger', ...
    'divdiff', ...
    'nintmva','nc.oi.exact','sdr','sdr.mva'}, ...
    'exact', 'deterministic');
% mixed open/closed limited load dependence: the Bruell-Balbo-Afshari
% effective-capacity MVA evaluates the product form itself, no expansion
add({'ncldmx'}, 'exact', 'deterministic');
add({'mm1','mmk','mxm1','mm1k','mg1','gm1','mapm1ps','pas'}, ...
    'exact', 'deterministic');
add({'mg1.prio','mg1.fb','mg1.srpt','mg1.psjf','mg1.setf','mg1.lrpt', ...
    'mm1.dps'}, 'exact', 'deterministic');
% MDD-rec: the normalising constant of a product form, summed EXACTLY over the
% reachable set by one memoised walk of the decision diagram that holds it. On a
% Petri net (solver_nc_spn_analyzer) and on a loss network (lossn_rec) alike.
add({'rec','lossn.rec','mdd.rec'}, 'exact', 'deterministic');
% state-space enumeration: every CTMC path but the perfect samplers
add({'ctmc','sync','flat','gpu','fd','uniformization'}, 'exact', 'deterministic');
% single-queue matrix-analytic results
add({'exact.mapmap1'}, 'exact', 'deterministic');
% JMT's exact analytical engine
add({'jmva','jmva.mva','jmva.recal','jmva.comom'}, 'exact', 'deterministic');
% loss networks
add({'lossn.exact'}, 'exact', 'deterministic');
% discrete-time (slotted) product form: the Bernoulli server of chapter 2 and
% the closed cycle of chapter 3 in Daduna (2001), both closed form
add({'dt.bernoulli1','dt.cycle','dt.cycleld'}, 'exact', 'deterministic');

% ---- exact, randomized ---------------------------------------------------
% coupling from the past samples the stationary law itself
add({'cftp','ctmc.cftp'}, 'exact', 'randomized');

% ---- approximate, deterministic ------------------------------------------
% approximate MVA and its variants
add({'amva','bs','aql','qsa','lin','gflin','egflin','dmlin','qd','qdlin','qdaql', ...
    'qli','fli','ab','schmidt','schmidt-ext','schmidtext','sqni','sum', ...
    'esum','cl','chandy-lakshmi','shadow','seidmann','linearizerms', ...
    'conway','rolia','zhou','suri','reiser.ms','chow','marie','sqd','mapqn', ...
    'lcp','pamb','pami','pamt','clust', ...
    'interp','highvar','balanced','tay','scat'}, 'approximate', 'deterministic');
% open-network decomposition and general-service closed forms
add({'qna','rqna','gig1','gigk','klb','kraemer','mg1k','mm1k.approx'}, ...
    'approximate', 'deterministic');
% asymptotic expansions and entropy/fixed-point methods over the
% normalizing constant
add({'le','ble','aghq','cub','kt','bkt','lekt','pana','panald','mem','mem.blocking','gm', ...
    'erlangfp','propfair','fpi','spm','ttl'}, ...
    'approximate', 'deterministic');
% balanced fairness aggregation is not closed under composition
add({'oi','balancedfairness','stationtime'}, 'approximate', 'deterministic');
% mean-field, diffusion and ODE methods
add({'fld','fluid','matrix','closing','statedep','softmin','pnorm','mfq', ...
    'rmf','tbi','diffusion','kp','minnormal','refined','dae'}, 'approximate', 'deterministic');
% phase-type network decomposition
add({'mam','dec','mna','ldqbd','qbd','qiu', ...
    'cdf','reneging','retrial','bgchain'}, 'approximate', 'deterministic');
% agent-based (RCAT) decomposition; 'exact' is the vestigial autocat alias and
% falls back to inap, so it is approximate like the rest of the family
add({'ag','inap','inapplus','inapinf'}, 'approximate', 'deterministic');
% layered and environment decomposition
add({'ln','ln.mva','layers','ln.dec','enhanced','ln.fluid','moment3','lqns','srvn', ...
    'srvn.ph','srvn.cs','flat','flat.cs','flat.ph','lqnsdefault','exactmva','srvn.exactmva','qns','env','env.blend', ...
    'blend','dec.avg'}, 'approximate', 'deterministic');
% approximate JMT analytical algorithms
add({'jmva.amva','jmva.chow','jmva.bs','jmva.aql','jmva.lin','jmva.dmlin'}, ...
    'approximate', 'deterministic');
% solver selection is a meta-method; the banner of the selected solver
% carries the real classification
add({'auto','tree','auto.tree','forest','cart'}, 'approximate', 'deterministic');

% ---- bound, deterministic ------------------------------------------------
% formal one-sided bounds; 'cub.upper', 'qrf.mem' and 'qrf.bas.mem' are spelt
% out because their tails 'cub' and 'mem' are NC method names
add({'ba','aba','bjb','mbjb','gb','pb','sb','mwba','pbh','pbk','bjbk','cbh', ...
    'ssd','sib','scb','ldbcmp','looping','bpt','bgt','qr','lr','qrf','cub.upper', ...
    'qrf.mem','harel', ...
    'auto.upper','auto.lower', ...
    'qrf.bas.mem','spnlp'}, 'bound', 'deterministic');

% ---- approximate, randomized ---------------------------------------------
% discrete-event simulation
add({'ssa','ldes','serial','para','parallel','nrm','jsim','replication', ...
    'jmt','lqsim','sim','uq'}, 'approximate', 'randomized');
% Monte Carlo and sampling-based evaluation of the normalizing constant
add({'mci','imci','amci','lhsmci','ls','is','sampling','lossn.mci'}, ...
    'approximate', 'randomized');
% Markov chain Monte Carlo on the regularized network (Chen-O'Cinneide)
add({'mcmc','nc.mcmc'}, 'approximate', 'randomized');
% the approximate sampler stops before coalescence
add({'cftp.approx'}, 'approximate', 'randomized');

% ---- per-solver default, for a method with no entry of its own -----------
% '#' prefixed so a solver name cannot be reached as a method name
add({'#ctmc'}, 'exact', 'deterministic');
add({'#ssa','#ldes','#jmt'}, 'approximate', 'randomized');
add({'#ba'}, 'bound', 'deterministic');
add({'#mva','#nc','#fld','#mam','#ln','#env','#lqns','#qns','#auto', ...
    '#uq'}, 'approximate', 'deterministic');
end
