function [runtime, analyzer] = runAnalyzer(self, options)
% RUNTIME = RUN()
% Run the solver

T0=tic;
if nargin<2
    options = self.getOptions;
end

QN = []; UN = [];
RN = []; TN = [];
CN = []; XN = [];
lG = NaN;

% The shared gate subsumes the coarse supports(self.model) test (through
% supportsModelMethod) and additionally reports which features are unsupported;
% a separate coarse check here would fire first and name nothing. It also
% validates options.method, as in the other solvers.
verboseGuard = self.runAnalyzerChecks(options); %#ok<NASGU> restores the caller verbosity on return

Solver.resetRandomGeneratorSeed(options.seed);

sn = getStruct(self); % doesn't need initial state
% The gate's own sentence for a caller running with enableChecks off: neither
% path can keep a self-looping job on its server (qns_immfeed_refusal).
immfeedReason = qns_immfeed_refusal(sn);
if ~isempty(immfeedReason)
    line_error(mfilename, immfeedReason);
end
line_ack('QNS', options.verbose);

line_debug(options, 'QNS: starting (method=%s, multiserver=%s)', options.method, options.config.multiserver);

method = options.method;

switch options.method
    case 'conway'
        options.config.multiserver = 'conway';
    case 'rolia'
        options.config.multiserver = 'rolia';
    case 'zhou'
        options.config.multiserver = 'zhou';
    case 'suri'
        options.config.multiserver = 'suri';
    case 'reiser'
        options.config.multiserver = 'reiser';
    case 'schmidt'
        options.config.multiserver = 'schmidt';
    case 'default'
        options.config.multiserver = 'rolia';
end

if self.model.hasProductFormSolution() || self.model.hasOpenClasses()
    % Use qnsolver directly for product-form networks or open networks
    % (QN2LQN does not support Source/Sink nodes for open networks)
    line_debug(options, 'QNS: product-form or open model, routing to qns_analyzer');
    [QN,UN,RN,TN,CN,XN,runtime,actualMethod] = solver_qns_analyzer(sn, options);
else
    line_debug(options, 'QNS: non-product-form closed model, converting to LQN and using SolverLQNS');
    lqnmodel=QN2LQN(self.model);
    lqn = lqnmodel.getStruct;
    tic;
    lqnsoptions = SolverLQNS.defaultOptions;
    lqnsoptions.verbose = false;
    actualMethod = options.method; % Track the actual method for LQNS path
    switch options.method
        case 'conway'
            lqnsoptions.config.multiserver = 'conway';
            actualMethod = 'conway';
        case 'rolia'
            lqnsoptions.config.multiserver = 'rolia';
            actualMethod = 'rolia';
        case 'zhou'
            lqnsoptions.config.multiserver = 'zhou';
            actualMethod = 'zhou';
        case 'suri'
            lqnsoptions.config.multiserver = 'suri';
            actualMethod = 'suri';
        case 'reiser'
            lqnsoptions.config.multiserver = 'reiser';
            actualMethod = 'reiser';
        case 'schmidt'
            lqnsoptions.config.multiserver = 'schmidt';
            actualMethod = 'schmidt';
        case 'default'
            actualMethod = 'rolia'; % Default is rolia for LQNS
    end
    AvgTable = SolverLQNS(lqnmodel,lqnsoptions).getAvgTable;
    runtime=toc;
    for r=1:sn.nclasses
        for i=1:sn.nstations
            t = lqn.ashift + r + (i-1)*sn.nclasses;
            QN(i,r) = AvgTable.QLen(t);
            % Convention crossing: lqns reports the processor utilization summed
            % over the host's servers, a Network station reports the per-server
            % fraction, and QN2LQN maps nservers to the host multiplicity.
            if ~isinf(sn.nservers(i))
                UN(i,r) = AvgTable.Util(t)/sn.nservers(i);
            else
                UN(i,r) = AvgTable.Util(t);
            end
            RN(i,r) = AvgTable.RespT(t);
            WN(i,r) = AvgTable.ResidT(t);
            TN(i,r) = AvgTable.Tput(t);
        end
    end
    XN=[];
    CN=[];
end

if nargout > 1
    analyzer = @(sn) solver_qns_analyzer(sn, options);
end

sn = self.getStruct;
T = getAvgTputHandles(self);
AN = sn_get_arvr_from_tput(sn, TN, T);

% Apply default(...) convention based on actual method used
if strcmp(method, 'default')
    method = ['default/' actualMethod];
else
    method = actualMethod;
end

self.setAvgResults(QN,UN,RN,TN,AN,[],CN,XN,runtime,method);

runtime = toc(T0);
end