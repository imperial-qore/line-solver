function [Q,U,R,T,C,X,lG,runtime,iter,method] = solver_nc_analyzer(sn, options)
% [Q,U,R,T,C,X,LG,RUNTIME,METHOD] = SOLVER_NC_ANALYZER(QN, OPTIONS)

% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.
Tstart = tic;

line_debug('NC analyzer starting: method=%s, nstations=%d, nclasses=%d, njobs=%s', options.method, sn.nstations, sn.nclasses, mat2str(sn.njobs));

% Order-independent (OI) closed network: explicit request, or auto-detected
% when the closed model consists solely of OI and delay stations. Solved
% exactly by the balanced-fairness normalizing constant (pfqn_ncoi) with the
% OI functional-server (pfqn_oi_fnc) identity for the mean queue lengths.
if nc_is_oi_model(sn) && any(strcmpi(options.method, {'default','exact'}))
    line_debug('NC analyzer routing to solver_nc_oi_analyzer (order-independent)');
    [Q,U,R,T,C,X,lG,runtime,iter,method] = solver_nc_oi_analyzer(sn, options);
    return
end

% IS specialized to OI/P&S stations -> pfqn_pas_is (reduces to pfqn_oi_is for an
% empty swap graph); see _kb/06-solver-catalog.md (NC section, analyzer routing)
if nc_is_pas_model(sn) && any(strcmpi(options.method, {'default','is','sampling'}))
    line_debug('NC analyzer routing to solver_nc_pas_is_analyzer (pass-and-swap IS)');
    [Q,U,R,T,C,X,lG,runtime,iter,method] = solver_nc_pas_is_analyzer(sn, options);
    return
end

% plain 'is' is the sample-an-ordering family for closed product-form networks;
% no open-class form. see _kb/06-solver-catalog.md (NC section, analyzer routing)
if strcmpi(options.method, 'is') && any(isinf(sn.njobs))
    line_error(mfilename, 'The ''is'' importance-sampling method requires a closed queueing network. Use ''sampling'' (pfqn_mci/pfqn_ls) for open or mixed models.');
end

% Maximum Entropy Method (Kouvatsos 1994): explicit request only. The
% method='default' path routes to the native normalizing-constant analyzer
% (solver_nc), as it did before MEM was introduced.
if strcmpi(options.method, 'mem')
    line_debug('NC analyzer routing to solver_nc_mem (Maximum Entropy)');
    [Q,U,R,T,C,X,iter,method] = solver_nc_mem(sn, options);
    lG = NaN;
    runtime = toc(Tstart);
    return
end

nservers = sn.nservers;
if max(nservers(nservers<Inf))>1 & any(isinf(sn.njobs)) & strcmpi(options.method,'exact') %#ok<AND2>
    line_error(mfilename,'NC solver cannot provide exact solutions for open or mixed queueing networks. Remove the ''exact'' option.');
end

% interpolate for non-integer closed populations
eta = abs(sn.njobs - floor(sn.njobs));
if any(eta>GlobalConstants.FineTol)
    line_debug('Fractional populations detected, using interpolation');
    sn_floor = sn; sn_floor.njobs = floor(sn.njobs);
    [Qf,Uf,Rf,Tf,Cf,Xf,lGf,~,iterf] = solver_nc(sn_floor, options);
    sn_ceil = sn; sn_ceil.njobs = ceil(sn.njobs);
    [Qc,Uc,Rc,Tc,Cc,Xc,lGc,~,iterc,method] = solver_nc(sn_ceil, options);
    Q = Qf +  eta .* (Qc-Qf);
    U = Uf +  eta .* (Uc-Uf);
    R = Rf +  eta .* (Rc-Rf);
    T = Tf +  eta .* (Tc-Tf);
    C = Cf +  eta .* (Cc-Cf);
    X = Xf +  eta .* (Xc-Xf);
    lG = lGf +  eta .* (lGf-lGc);
    iter = iterc + iterf;
    line_debug('NC interpolation complete: used %d + %d iterations', iterf, iterc);
else % if integers or open model
    if any(isinf(sn.njobs))
        line_debug('Open/mixed model detected, calling solver_nc');
    else
        line_debug('Using exact integer populations, calling solver_nc');
    end
    [Q,U,R,T,C,X,lG,~,iter,method] = solver_nc(sn, options);
end
runtime = toc(Tstart);
end
