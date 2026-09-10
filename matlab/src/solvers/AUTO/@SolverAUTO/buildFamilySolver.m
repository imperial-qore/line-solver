function solver = buildFamilySolver(self, family, model, options) %#ok<INUSL>
% SOLVER = BUILDFAMILYSOLVER(FAMILY, MODEL, OPTIONS)
%
% Instantiate the solver of a method family. OPTIONS.METHOD already holds the
% submethod resolved by RESOLVEMETHODTOKEN, so a pinned method reaches the
% family that runs it.

switch family
    case 'mam'
        solver = SolverMAM(model, options);
    case 'ag'
        solver = SolverAG(model, options);
    case 'mva'
        solver = SolverMVA(model, options);
    case 'nc'
        solver = SolverNC(model, options);
    case 'fluid'
        solver = SolverFluid(model, options);
    case 'jmt'
        solver = SolverJMT(model, options);
    case 'ssa'
        solver = SolverSSA(model, options);
    case 'ctmc'
        solver = SolverCTMC(model, options);
    case 'ldes'
        solver = SolverLDES(model, options);
    case 'ba'
        solver = SolverBA(model, options);
    case 'qns'
        solver = SolverQNS(model, options);
    case 'lqns'
        solver = SolverLQNS(model, options);
    case 'ln'
        [factory, outer] = innerFactory(options);
        solver = SolverLN(model, factory, outer);
    case 'env'
        [factory, outer] = innerFactory(options);
        solver = SolverENV(model, factory, outer);
    case 'uq'
        inner = options;
        inner.method = 'default';
        solver = SolverUQ(model, @(m) LINE(m, inner), options);
    otherwise
        line_error(mfilename, sprintf('Unknown method family ''%s''.', family));
end
end

function [factory, outer] = innerFactory(options)
% [FACTORY, OUTER] = INNERFACTORY(OPTIONS)
%
% A composite family may be qualified by the family that solves its inner
% models, as in 'ln.mva' or 'env.fluid'. Anything else is an intent, and the
% inner models are then solved by LINE itself, which selects per submodel.

outer = options;
inner = options;
fam = SolverAUTO.familyAlias(options.method);
if ~isempty(fam) && ~strcmp(fam,'ln') && ~strcmp(fam,'env')
    outer.method = 'default';
    inner.method = 'default';
    inner.verbose = 0;
    switch fam
        case 'mva'
            factory = @(m) SolverMVA(m, inner);
        case 'nc'
            factory = @(m) SolverNC(m, inner);
        case 'mam'
            factory = @(m) SolverMAM(m, inner);
        case 'fluid'
            factory = @(m) SolverFluid(m, inner);
        case 'ctmc'
            factory = @(m) SolverCTMC(m, inner);
        case 'ssa'
            factory = @(m) SolverSSA(m, inner);
        case 'ldes'
            factory = @(m) SolverLDES(m, inner);
        case 'jmt'
            factory = @(m) SolverJMT(m, inner);
        otherwise
            factory = @(m) LINE(m, inner);
    end
else
    inner.verbose = 0;
    factory = @(m) LINE(m, inner);
end
end
