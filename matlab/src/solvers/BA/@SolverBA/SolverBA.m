classdef SolverBA < NetworkSolver
    % Bound Analysis solver for closed queueing networks.
    %
    % SolverBA is a dedicated home for asymptotic and hierarchical throughput/
    % queue-length BOUNDS. Unlike SolverMVA (point estimates), each method
    % returns an optimistic or pessimistic bound; the .upper/.lower pair for a
    % family brackets the exact solution. Bounds need only demands (visits x
    % service time) and populations -- no service-time distributions -- so the
    % featset is deliberately narrow (closed product-form-parameterized models).
    %
    % Method families:
    %   Noniterative (SolverBA-native, solver_ba_analyzer):
    %     aba.*  Asymptotic Bound Analysis (Denning-Buzen)
    %     bjb.*  Balanced Job Bounds (Zahorjan et al.)
    %     pb.*   Proportional Bounds (Eager-Sevcik)
    %     gb.*   Geometric Bounds (Casale-Muntz-Serazzi)
    %     sb.*   Simple bounds from power sums (Harel-Namn-Sturm); NOT a
    %            geometric variant despite sitting beside gb. Rejects delay
    %            stations, as does lr; gb accepts them.
    %     mwba.* Multiclass worst-case balanced bound
    %   Hierarchical / iterative (SolverBA-native, level-parameterized):
    %     pbh    Performance Bound Hierarchy (Eager-Sevcik 1983), level option
    %     sib    Successively Improving Bounds (Srinivasan 1985), level option
    %     cbh    Convolutional Bound Hierarchies (Dowdy et al. 1984)
    %     bjbk / pbk  iterative BJB(k)/PB(k) with delay (Casale et al. 2008)
    %     cub    Composite Upper Bound, multiclass upper-only (Kerola 1986)
    %     mbjb   multiclass Balanced Job Bounds lower (Kerola eq. 10; seeds cub)
    %     ssd    multiServer Disaggregation bounds (Suri-Dallery 1986)
    %     ldbcmp LD-BCMP closed-open equivalence bound (Anselmi-Cremonesi 2008)
    %   LP-based reduction bounds (QRF library, solver_ba_qrf_analyzer):
    %     qr / qrf.mmi        Quadratic Reduction Framework (single-class PH)
    %     lr / qrf.mmi.linear Linear Reduction variant
    %     qrf.mem, qrf.mmi.ld, qrf.bas.*, qrf.rsrd  further QRF sub-methods
    %
    % Copyright (c) 2012-2026, Imperial College London
    % All rights reserved.

    methods
        function self = SolverBA(model, varargin)
            % SOLVERBA Create a bound-analysis solver instance
            self@NetworkSolver(model, mfilename);
            self.setOptions(Solver.parseOptions(varargin, SolverBA.defaultOptions));
            self.setLang();
        end

        function sn = getStruct(self)
            % GETSTRUCT Get model data structure for analysis
            sn = self.model.getStruct(false);
        end

        [runtime, analyzer] = runAnalyzer(self, options);

        function bounds = getBounds(self)
            % GETBOUNDS Return the {lower,upper} bracket for the requested
            % bound family as a struct with fields Xlower/Xupper (throughput)
            % and Qlower/Qupper (queue lengths). The current method's family
            % prefix (before the first '.') selects the family; a hierarchical
            % method uses its own level option.
            % One-sided families (cub upper-only, mbjb/ldbcmp lower-only)
            % return NaN on the missing side.
            fam = self.options.method;
            dot = strfind(fam,'.');
            if ~isempty(dot)
                fam = fam(1:dot(1)-1);
            end
            valid = self.listValidMethods();
            Tl = NaN; Ql = NaN; Tu = NaN; Qu = NaN;
            if any(strcmp([fam,'.lower'], valid))
                [Ql,~,~,Tl] = self.solveSide([fam,'.lower']);
            end
            if any(strcmp([fam,'.upper'], valid))
                [Qu,~,~,Tu] = self.solveSide([fam,'.upper']);
            end
            bounds = struct('Tlower',Tl,'Tupper',Tu,'Qlower',Ql,'Qupper',Qu);
        end

        function [Q,U,R,T] = solveSide(self, method)
            % SOLVESIDE Re-run the solver for one side of the bracket.
            % The re-run instance inherits the CALLER'S FULL OPTION SET (level,
            % verbose, tol, ...) and only overrides the method. Constructing it
            % with just 'method' would silently reset options.level to its
            % default of 2, so a hierarchical family (pbh/cbh/pbk/bjbk/sib)
            % reached through getBounds would never tighten as level is raised.
            s = SolverBA(self.model);
            opts = self.options;
            opts.method = method;
            s.setOptions(opts);
            [Q,U,R,T] = s.getAvg();
        end

        function BoundsTable = getBoundsTable(self, keepDisabled)
            % GETBOUNDSTABLE Table of the {lower,upper} bracket per station and
            % class, in the layout of GETAVGTABLE. Columns:
            %   Station, JobClass, Qlower, Qupper, Tlower, Tupper
            % One-sided families (cub upper-only, mbjb/ldbcmp lower-only) carry
            % NaN on the missing side; NaN is preserved, never replaced by zero.
            if nargin < 2
                keepDisabled = false;
            end
            b = self.getBounds();
            sn = self.getStruct();
            M = sn.nstations;
            K = sn.nclasses;
            Ql = SolverBA.expandBound(b.Qlower, M, K);
            Qu = SolverBA.expandBound(b.Qupper, M, K);
            Tl = SolverBA.expandBound(b.Tlower, M, K);
            Tu = SolverBA.expandBound(b.Tupper, M, K);
            [Qlval, Quval, Tlval, Tuval] = deal([]);
            JobClass = {};
            Station = {};
            for ist = 1:M
                for k = 1:K
                    vals = [Ql(ist,k), Qu(ist,k), Tl(ist,k), Tu(ist,k)];
                    finite = vals(~isnan(vals));
                    % Mirror getAvgTable's drop of disabled station-class pairs,
                    % but NaN-safe: a row is kept when any value that is present
                    % is nonzero, so an all-NaN side never removes the row.
                    if keepDisabled || isempty(finite) || any(finite ~= 0)
                        JobClass{end+1,1} = sn.classnames{k}; %#ok<AGROW>
                        Station{end+1,1} = sn.nodenames{sn.stationToNode(ist)}; %#ok<AGROW>
                        Qlval(end+1) = Ql(ist,k); %#ok<AGROW>
                        Quval(end+1) = Qu(ist,k); %#ok<AGROW>
                        Tlval(end+1) = Tl(ist,k); %#ok<AGROW>
                        Tuval(end+1) = Tu(ist,k); %#ok<AGROW>
                    end
                end
            end
            Station = label(Station);
            JobClass = label(JobClass);
            Qlower = Qlval(:); Qupper = Quval(:);
            Tlower = Tlval(:); Tupper = Tuval(:);
            BoundsTable = Table(Station, JobClass, Qlower, Qupper, Tlower, Tupper);
            BoundsTable = IndexedTable(BoundsTable);
        end

        function [allMethods] = listValidMethods(self)
            % LISTVALIDMETHODS Valid bound methods for the current model.
            % Hierarchical methods (pbh/sib/cbh/bjbk/pbk/cub/ssd/ldbcmp) are
            % appended here as each backing pfqn_* algorithm is integrated.
            allMethods = { ...
                'default', ...
                'aba.upper','aba.lower', ...
                'bjb.upper','bjb.lower', ...
                'pb.upper','pb.lower', ...
                'gb.upper','gb.lower', ...
                'sb.upper','sb.lower', ...
                'mwba.upper','mwba.lower', ...
                'pbh.upper','pbh.lower', ...
                'pbk.upper','pbk.lower', ...
                'bjbk.upper','bjbk.lower', ...
                'cbh.upper','cbh.lower', ...
                'ssd.upper','ssd.lower', ...
                'cub.upper','mbjb.lower', ...
                'sib.upper','sib.lower', ...
                'ldbcmp.lower', ...
                'qr','lr','lr.upper','lr.lower', ...
                'qrf.mmi','qrf.mem','qrf.mmi.ld','qrf.mmi.linear', ...
                'qrf.bas.mmi','qrf.bas.mem','qrf.bas','qrf.rsrd'};
        end
    end

    methods (Static)
        function Y = expandBound(X, M, K)
            % EXPANDBOUND Normalize a bracket side to an (M,K) matrix. A
            % one-sided family leaves its missing side as the scalar NaN
            % returned by getBounds; expand it so the table keeps full shape
            % with NaN entries rather than erroring or dropping the column.
            if isempty(X)
                Y = NaN(M,K);
            elseif isscalar(X)
                Y = repmat(X, M, K);
            else
                Y = X;
            end
        end

        function libs = getLibrariesUsed(sn, options) %#ok<INUSL>
            % GETLIBRARIESUSED External libraries used by SolverBA. The QRF
            % reduction bounds solve an LP via the Optimization Toolbox.
            libs = {};
            if ~isempty(options) && isfield(options,'method') && ...
                    (startsWith(options.method,'qrf') || any(strcmp(options.method,{'qr','lr'})))
                libs{end+1} = 'MATLAB Optimization Toolbox (fmincon)';
            end
        end

        function options = defaultOptions
            % OPTIONS = DEFAULTOPTIONS()
            options = SolverOptions('MVA');
            options.method = 'default';
            % Hierarchy level for pbh/cbh (MVA steps / exact-convolved servers)
            % and iteration count k for pbk/bjbk. Default 2.
            options.level = 2;
        end

        function [bool, featSupported] = supports(model)
            % SUPPORTS Whether SolverBA can bound the given model
            featUsed = model.getUsedLangFeatures();
            featSupported = SolverFeatureSet();
            featSupported.setTrue({ ...
                'ClassSwitch','DelayStation','Queue', ...
                'Sink','Source','Router', ...
                'ClosedClass','OpenClass', ...
                'SchedStrategy_INF','SchedStrategy_PS', ...
                'SchedStrategy_FCFS','SchedStrategy_LCFSPR', ...
                'RoutingStrategy_PROB','RoutingStrategy_RAND', ...
                'ClosedClass_multiclass'});
            bool = SolverFeatureSet.supports(featSupported, featUsed);
        end
    end
end
