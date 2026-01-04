function [runtime, tranSysState, tranSync] = runAnalyzer(self, options)
% [RUNTIME, TRANSYSSTATE] = RUNANALYZER()
%
% Run the DES solver on the queueing network model.
% For LayeredNetwork models, use SolverLDES from line-apps.

T0 = tic;
if nargin < 2
    options = self.getOptions;
end
self.runAnalyzerChecks(options);
Solver.resetRandomGeneratorSeed(options.seed);

line_debug('DES solver starting: lang=%s, samples=%d, seed=%d', options.lang, options.samples, options.seed);

tranSysState = [];
tranSync = [];

% Check if confidence intervals are requested
[confintEnabled, confintLevel] = Solver.parseConfInt(options.confint);

switch options.lang
    case 'java'
        line_debug('Default method: using Java DES discrete-event simulation\n');
        line_debug('Using Java DES backend');
        % Convert model to Java
        jmodel = LINE2JLINE(self.model);

        % Create Java solver with DES-specific options via JLINE wrapper
        jsolver = JLINE.SolverDES(jmodel, options);

        % Regular Network analysis
        [runtime] = runNetworkAnalyzer(self, jmodel, jsolver, confintEnabled, options);
    otherwise
        line_error(mfilename, 'SolverDES currently only supports Java backend. Use options.lang = ''java''.');
end

runtime = toc(T0);
end

function [runtime] = runNetworkAnalyzer(self, jmodel, jsolver, confintEnabled, options)
% RUNNETWORKANALYZER Run DES analysis for regular Network models

M = jmodel.getNumberOfStations;
R = jmodel.getNumberOfClasses;

% Check if this is transient analysis
isTransient = isfield(options, 'timespan') && length(options.timespan) >= 2 && ...
              isfinite(options.timespan(2));

if isTransient
    % Run transient analysis directly
    runTransientAnalyzer(self, jsolver, M, R, options);
    runtime = 0;  % Runtime tracked inside transient analyzer
    return;
end

% Steady-state analysis
[QN, UN, RN, WN, AN, TN] = JLINE.arrayListToResults(jsolver.getAvgTable(true));

CN = JLINE.from_jline_matrix(jsolver.getAvgSysRespT());
XN = JLINE.from_jline_matrix(jsolver.getAvgSysTput());
runtime = jsolver.result.runtime;
QN = reshape(QN', R, M)';
UN = reshape(UN', R, M)';
RN = reshape(RN', R, M)';
TN = reshape(TN', R, M)';
WN = reshape(WN', R, M)';
AN = reshape(AN', R, M)';

% Extract FCR metrics and append to result matrices
sn = self.model.getStruct;
F = sn.nregions;
if F > 0
    result = jsolver.result;
    isValidMatrix = @(x) ~isempty(x) && isa(x, 'jline.util.matrix.Matrix');
    if isValidMatrix(result.QNfcr)
        QNfcr = JLINE.from_jline_matrix(result.QNfcr);
        UNfcr = JLINE.from_jline_matrix(result.UNfcr);
        RNfcr = JLINE.from_jline_matrix(result.RNfcr);
        TNfcr = JLINE.from_jline_matrix(result.TNfcr);
        WNfcr = JLINE.from_jline_matrix(result.WNfcr);
        % ANfcr is NaN for FCR (not applicable)
        ANfcr = NaN(F, R);
        % Append FCR rows to station metrics
        QN = [QN; QNfcr];
        UN = [UN; UNfcr];
        RN = [RN; RNfcr];
        TN = [TN; TNfcr];
        WN = [WN; WNfcr];
        AN = [AN; ANfcr];
    end
end

% Print samples like SSA does
if options.verbose
    line_printf('DES samples: %8d\n', options.samples);
end

self.setAvgResults(QN, UN, RN, TN, AN, WN, CN, XN, runtime, options.method, options.samples);

% Extract confidence intervals from Java solver results
if confintEnabled
    result = jsolver.result;
    % Helper to check for non-null Java objects
    isValidMatrix = @(x) ~isempty(x) && isa(x, 'jline.util.matrix.Matrix');
    if isValidMatrix(result.QNCI)
        QNCI = JLINE.from_jline_matrix(result.QNCI);
        QNCI = reshape(QNCI', R, M)';
    else
        QNCI = [];
    end
    if isValidMatrix(result.UNCI)
        UNCI = JLINE.from_jline_matrix(result.UNCI);
        UNCI = reshape(UNCI', R, M)';
    else
        UNCI = [];
    end
    if isValidMatrix(result.RNCI)
        RNCI = JLINE.from_jline_matrix(result.RNCI);
        RNCI = reshape(RNCI', R, M)';
    else
        RNCI = [];
    end
    if isValidMatrix(result.TNCI)
        TNCI = JLINE.from_jline_matrix(result.TNCI);
        TNCI = reshape(TNCI', R, M)';
    else
        TNCI = [];
    end
    if isValidMatrix(result.ANCI)
        ANCI = JLINE.from_jline_matrix(result.ANCI);
        ANCI = reshape(ANCI', R, M)';
    else
        ANCI = [];
    end
    if isValidMatrix(result.WNCI)
        WNCI = JLINE.from_jline_matrix(result.WNCI);
        WNCI = reshape(WNCI', R, M)';
    else
        WNCI = [];
    end
    % Store CI results
    self.setAvgResultsCI(QNCI, UNCI, RNCI, TNCI, ANCI, WNCI, [], []);
end

end

function runTransientAnalyzer(self, jsolver, Mjava, R, options)
% RUNTRANSIENTANALYZER Run transient analysis for DES

T0 = tic;

% Get MATLAB model structure for correct station count
sn = self.model.getStruct;
M = sn.nstations;  % Use MATLAB station count for cell arrays

% Run transient analysis on Java side
jsolver.getTranAvg();
result = jsolver.result;
runtime = toc(T0);

% Check if transient results are available (result.QNt is a Java 2D array)
if ~isempty(result.QNt) && ~isempty(result.t)
    % Get number of time points
    Tmax = result.t.length();
    if Tmax > 0
        % Convert Java transient results to MATLAB format
        % result.QNt(ist,r) is a Matrix(numTimePoints, 2) with columns [value, time]
        % Use MATLAB station count for cell arrays so getTranAvg can index correctly
        QNt = cell(M, R);
        UNt = cell(M, R);
        RNt = cell(M, R);  % Not computed by DES, will be empty
        TNt = cell(M, R);
        CNt = cell(1, R);
        XNt = cell(1, R);

        % Initialize all cells to NaN
        for ist = 1:M
            for r = 1:R
                QNt{ist, r} = NaN;
                UNt{ist, r} = NaN;
                RNt{ist, r} = NaN;
                TNt{ist, r} = NaN;
            end
        end

        % Extract data from Java results (indexed by Java station count)
        for ist = 1:Mjava
            for r = 1:R
                % Queue length transient
                try
                    jmatrix = result.QNt(ist, r);
                    if ~isempty(jmatrix) && jmatrix.getNumRows() > 0
                        nRows = jmatrix.getNumRows();
                        matData = NaN(nRows, 2);
                        for p = 1:nRows
                            matData(p, 1) = jmatrix.get(p-1, 0);  % value
                            matData(p, 2) = jmatrix.get(p-1, 1);  % time
                        end
                        QNt{ist, r} = matData;
                    else
                        QNt{ist, r} = NaN;
                    end
                catch
                    QNt{ist, r} = NaN;
                end

                % Utilization transient
                try
                    jmatrix = result.UNt(ist, r);
                    if ~isempty(jmatrix) && jmatrix.getNumRows() > 0
                        nRows = jmatrix.getNumRows();
                        matData = NaN(nRows, 2);
                        for p = 1:nRows
                            matData(p, 1) = jmatrix.get(p-1, 0);  % value
                            matData(p, 2) = jmatrix.get(p-1, 1);  % time
                        end
                        UNt{ist, r} = matData;
                    else
                        UNt{ist, r} = NaN;
                    end
                catch
                    UNt{ist, r} = NaN;
                end

                % Throughput transient
                try
                    jmatrix = result.TNt(ist, r);
                    if ~isempty(jmatrix) && jmatrix.getNumRows() > 0
                        nRows = jmatrix.getNumRows();
                        matData = NaN(nRows, 2);
                        for p = 1:nRows
                            matData(p, 1) = jmatrix.get(p-1, 0);  % value
                            matData(p, 2) = jmatrix.get(p-1, 1);  % time
                        end
                        TNt{ist, r} = matData;
                    else
                        TNt{ist, r} = NaN;
                    end
                catch
                    TNt{ist, r} = NaN;
                end

                % Response time transient not computed by DES
                RNt{ist, r} = NaN;
            end
        end

        % CNt and XNt not computed for transient
        for r = 1:R
            CNt{1, r} = NaN;
            XNt{1, r} = NaN;
        end

        % Store transient results
        self.setTranAvgResults(QNt, UNt, RNt, TNt, CNt, XNt, runtime);

        % Also compute and store steady-state estimates from final transient values
        QN = NaN(M, R);
        UN = NaN(M, R);
        RN = NaN(M, R);
        TN = NaN(M, R);
        WN = NaN(M, R);
        AN = NaN(M, R);
        CN = NaN(1, R);
        XN = NaN(1, R);

        for ist = 1:M
            for r = 1:R
                if ~isscalar(QNt{ist, r}) && ~isnan(QNt{ist, r}(end, 1))
                    QN(ist, r) = QNt{ist, r}(end, 1);
                end
                if ~isscalar(UNt{ist, r}) && ~isnan(UNt{ist, r}(end, 1))
                    UN(ist, r) = UNt{ist, r}(end, 1);
                end
                if ~isscalar(TNt{ist, r}) && ~isnan(TNt{ist, r}(end, 1))
                    TN(ist, r) = TNt{ist, r}(end, 1);
                end
            end
        end

        self.setAvgResults(QN, UN, RN, TN, AN, WN, CN, XN, runtime, options.method, 0);
    end
end

% Print verbose output if requested
if options.verbose
    line_printf('DES transient analysis complete, timespan = [%g, %g]\n', ...
        options.timespan(1), options.timespan(2));
end
end
