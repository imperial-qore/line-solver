classdef ParamEstimator < handle
    % Abstract class for service demand estimators
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.
    
    properties (Access = public)
        options; % Data structure with engine options
        model; % Model to be analyzed
        samples; % Input dataset
        samplesAggr; % Input dataset
    end
    
    methods (Hidden)
        %Constructor
        function self = ParamEstimator(model, options)
            % SELF = SOLVER(MODEL, NAME, OPTIONS)
            if ~exist('options','var')
                options = self.defaultOptions;
            end
            self.model = model;
            self.options = options;
            self.samples = cell(model.getNumberOfNodes, model.getNumberOfClasses);
            self.samplesAggr = cell(model.getNumberOfNodes, 1);
        end
    end
    
    methods
        function self = addSamples(self, sampleData)
            i = self.model.getNodeIndex(sampleData.node);
            r = 0;
            if sampleData.isAggregate
                %jobclass = NaN;
                self.samplesAggr{i}{end+1} = sampleData;
            else
                r = self.model.getClassIndex(sampleData.class);
                % store
                self.samples{i,r}{end+1} = sampleData;
            end
        end
        
        function data = getData(self)
            data = self.samples;
        end
        
        function data = getDataAggr(self)
            data = self.samplesAggr;
        end
        
        % look into the data available for that node and class for the
        % first ArvR dataset
        function data = getArvR(self, node, jobclass)
            data = [];
            i = self.model.getNodeIndex(node);
            r = self.model.getClassIndex(jobclass);
            nodeData = self.samples{i,r};
            for d=1:length(nodeData)
                if nodeData{d}.type == MetricType.ArvR
                    data = nodeData{d};
                    return
                end
            end
        end
        
        % look into the data available for that node and class for the
        % first ArvR dataset
        function data = getUtil(self, node, jobclass)
            data = [];
            i = self.model.getNodeIndex(node);
            r = self.model.getClassIndex(jobclass);
            nodeData = self.samples{i,r};
            for d=1:length(nodeData)
                if nodeData{d}.type == MetricType.Util
                    data = nodeData{d};
                    return
                end
            end
        end

        % look into the data available for that node and class for the
        % first ArvR dataset
        function data = getRespT(self, node, jobclass)
            data = [];
            i = self.model.getNodeIndex(node);
            r = self.model.getClassIndex(jobclass);
            nodeData = self.samples{i,r};
            for d=1:length(nodeData)
                if nodeData{d}.type == MetricType.RespT
                    data = nodeData{d};
                    return
                end
            end
        end
        
        % look into the data available for that node and class for the
        % first ArvR dataset
        function data = getAggrUtil(self, node)
            data = [];
            i = self.model.getNodeIndex(node);
            nodeAggrData = self.samplesAggr{i};
            for d=1:length(nodeAggrData)
                if nodeAggrData{d}.type == MetricType.Util
                    data = nodeAggrData{d};
                    return
                end
            end
        end
        
        % look into the data available for that node and class for the
        % first QLen dataset
        function data = getQLen(self, node, jobclass, ev)
            data = [];
            i = self.model.getNodeIndex(node);
            r = self.model.getClassIndex(jobclass);
            nodeData = self.samples{i,r};
            if ~exist('ev', 'var')
                for d=1:length(nodeData)
                    if nodeData{d}.type == MetricType.QLen
                        data{end+1} = nodeData{d};
                    end
                end
                if length(data) == 1
                    data = data{1};
                end
            else
                for d=1:length(nodeData)
                    % e.g., arrival queue-length
                    if nodeData{d}.type == MetricType.QLen && (nodeData{d}.cond.node == ev.node) && (nodeData{d}.cond.class == ev.class) && (nodeData{d}.cond.event == ev.event)
                        data = nodeData{d};
                        return
                    end
                end
            end
        end
        
        % look into the data available for that node and class for the
        % first QLen dataset
        function data = getAggrQLen(self, node, ev)
            data = [];
            i = self.model.getNodeIndex(node);
            nodeData = self.samplesAggr{i};
            if ~exist('ev', 'var')
                for d=1:length(nodeData)
                    if nodeData{d}.type == MetricType.QLen
                        data = nodeData{d};
                        return
                    end
                end
            else
                for d=1:length(nodeData)
                    % e.g., arrival queue-length
                    if nodeData{d}.type == MetricType.QLen && (nodeData{d}.cond.node == ev.node) && (nodeData{d}.cond.class == ev.class) && (nodeData{d}.cond.event == ev.event)
                        data = nodeData{d};
                        return
                    end
                end
            end
        end

        % look into the data available for that node and class for the
        % first Tput dataset
        function data = getTput(self, node, jobclass)
            data = [];
            i = self.model.getNodeIndex(node);
            r = self.model.getClassIndex(jobclass);
            nodeData = self.samples{i,r};
            for d=1:length(nodeData)
                if nodeData{d}.type == MetricType.Tput
                    data = nodeData{d};
                    return
                end
            end
        end

        function method = autoMethod(self)
            % AUTOMETHOD Automatically select the best estimation method
            % based on the available sampled metrics.
            sn = self.model.getStruct;
            hasArvR = false; hasRespT = false; hasUtil = false;
            hasQLen = false; hasTput = false; hasTrace = false;
            hasAggrUtil = false; hasAggrQLen = false;

            % scan per-class metrics
            for i = 1:size(self.samples, 1)
                for r = 1:size(self.samples, 2)
                    nodeData = self.samples{i,r};
                    for d = 1:length(nodeData)
                        sm = nodeData{d};
                        if sm.type == MetricType.ArvR, hasArvR = true; end
                        if sm.type == MetricType.RespT, hasRespT = true; end
                        if sm.type == MetricType.Util, hasUtil = true; end
                        if sm.type == MetricType.QLen, hasQLen = true; end
                        if sm.type == MetricType.Tput, hasTput = true; end
                        if sm.isTrace(), hasTrace = true; end
                    end
                end
            end

            % scan aggregate metrics
            for i = 1:size(self.samplesAggr, 1)
                nodeAggrData = self.samplesAggr{i};
                for d = 1:length(nodeAggrData)
                    sm = nodeAggrData{d};
                    if sm.type == MetricType.Util, hasAggrUtil = true; end
                    if sm.type == MetricType.QLen, hasAggrQLen = true; end
                end
            end

            % select method based on available data
            if hasTrace && hasRespT && hasAggrQLen
                method = 'erps';
            elseif hasTrace && hasRespT && hasArvR
                method = 'mlps';
            elseif hasArvR && hasRespT && (hasUtil || hasAggrUtil)
                method = 'ubo';
            elseif hasArvR && (hasUtil || hasAggrUtil)
                method = 'ubr';
            elseif hasQLen
                method = 'qmle';
            else
                error('Insufficient data to automatically select an estimation method. Please set options.method manually.');
            end

            self.options.method = method;
        end

        interpolate(self);
        estVal = estimateAt(self, nodes);
        estVal = estimator_ubo(self, nodes);
        estVal = estimator_ubr(self, nodes);
        estVal = estimator_erps(self, nodes);
        estVal = estimator_ekf(self, nodes);
        estVal = estimator_mcmc(self, nodes);
        estVal = estimator_mle(self, nodes);
        estVal = estimator_rnn(self, nodes);
        estVal = estimator_mlps(self, nodes);
        estVal = estimator_fmlps(self, nodes);
        estVal = estimator_qmle(self, nodes);
        [eqModel, eqNode] = buildClosedEquivalentForPS(self, node);
        estVal = estimator_gibbs(self, nodes);
    end

    methods (Static)
        function options = defaultOptions()
            % OPTIONS = DEFAULTOPTIONS()
            % Return default options
            options = struct();
            options.verbose = 1;
            options.method = 'ubr';
            options.variant = 'default';
            options.iter_max = 1000;
            options.tol = 1e-3;
            options.solver = @SolverAuto;
            options.openPopulation = 100;
        end

        function desc = getRequiredMetrics(method)
            % GETREQUIREDMETRICS Return description of metrics required
            % by each estimation method.
            switch method
                case 'ubr'
                    desc = 'ArvR (per-class) + Util (per-class or aggregate)';
                case 'ubo'
                    desc = 'ArvR (per-class) + RespT (per-class) + Util (aggregate)';
                case 'erps'
                    desc = 'RespT (per-class) + QLen (aggregate, conditional on class arrivals). PS stations only.';
                case 'ekf'
                    desc = 'RespT (per-class) + Util (aggregate). Sequential/recursive estimation.';
                case 'mcmc'
                    desc = 'QLen (aggregate). Gibbs sampling with MCMC. Open/mixed via closed equivalence.';
                case 'mle'
                    desc = 'ArvR (per-class) + RespT (per-class) + Util (aggregate)';
                case 'rnn'
                    desc = 'QLen (per-class, trace format). Transient queue-length traces.';
                case 'mlps'
                    desc = 'ArvR (per-class, trace) + RespT (per-class, trace). PS stations only. Open/mixed via closed equivalence.';
                case 'fmlps'
                    desc = 'ArvR (per-class, trace) + RespT (per-class, trace). PS stations only. Open/mixed via closed equivalence.';
                case 'qmle'
                    desc = 'QLen (per-class). Open/mixed via closed equivalence (Z_r = N_r / lambda_r).';
                case 'gibbs'
                    desc = 'ArvR (per-class, trace) + RespT (per-class, trace) + Tput (per-class). Gibbs sampling.';
                otherwise
                    desc = sprintf('Unknown method: %s', method);
            end
        end
    end
end