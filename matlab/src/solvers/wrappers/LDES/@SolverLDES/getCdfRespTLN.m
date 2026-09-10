function RD = getCdfRespTLN(self)
% RD = GETCDFRESPTLN()
%
% Empirical response time distribution of every ENTRY of a LayeredNetwork,
% measured on the simulated sample path.
%
% The simulated counterpart of @SolverLN/getCdfRespT: where the moment3 pass
% fits an APH to three moments and convolves, this is the ecdf of the response
% times the run observed, so its tail is measured rather than extrapolated. The
% engine times each request from the instant the entry acquires a thread to the
% instant it replies -- the interval RLN averages -- so the mean of this law
% reproduces that row.
%
% RD is an (nentries x 1) cell of [F(t), t] matrices, in the entry-local index
% space (lsn.eshift + i), the column order every getCdfRespT in LINE follows. An
% entry the run observed nothing at is left EMPTY rather than filled with a
% guess.
%
% LayeredNetwork LDES runs on the Java ensemble backend (self.obj), not through
% the JSON CLI, so the samples are read straight off that object.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isa(self.model, 'LayeredNetwork')
    line_error(mfilename, ['getCdfRespTLN requires a LayeredNetwork; this solver holds a ' ...
        'Network, whose response times are indexed by (station, class). Use getCdfRespT.']);
end

lsn = self.model.getStruct();
RD = cell(lsn.nentries, 1);
if GlobalConstants.DummyMode
    return
end

T0 = tic;
if isempty(self.obj)
    self.setLang();
end
% The samples ride on the LN result, so the run has to have happened.
self.obj.getAvg();
jcdf = self.obj.getCdfRespTLN();

for e = 1:lsn.nentries
    cell_e = jcdf.get(e-1); % java.util.List is 0-based
    if isempty(cell_e)
        continue
    end
    RD{e} = JLINE.from_jline_matrix(cell_e);
end

runtime = toc(T0);
self.setDistribResults(RD, runtime);
end
