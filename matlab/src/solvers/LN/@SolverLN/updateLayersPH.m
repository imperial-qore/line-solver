function updateLayersPH(self, it) %#ok<INUSD>
% UPDATELAYERSPH(IT) Push the composed laws into the layers of a PH encoding
%
% A layer of these methods carries no routing that depends on the iterate: the
% number of calls a caller makes is folded into its service law rather than into
% a visit ratio, so only two laws move per (server, class) -- the phase-type
% service law at the server and the mean of the surrogate delay at the client.
%
% Serves 'srvn.ph' and 'flat.ph' alike. Everything the two share is written
% against ph.layer{idx}, whose station indices are absolute, and the one place
% they differ is phDelayMean: a surrogate delay stands for whatever of a
% thread's cycle the layer does NOT hold, and under the squashed layering it
% holds nearly all of it.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

% updateMetricsPH has already recomposed the entry laws from this iterate; a
% second pass here would only refit the same leaves
lqn = self.lqn;
ph = self.ph;

for idx = 1:(lqn.nhosts + lqn.ntasks)
    if isnan(self.idxhash(idx)) || isempty(ph.layer{idx})
        continue
    end
    L = ph.layer{idx};
    model = self.ensemble{self.idxhash(idx)};
    clientDelay = model.nodes{1};

    for c = L.callers
        k = L.classOfCaller(c);
        cls = model.classes{k};
        [alpha, T] = phServiceLaw(self, lqn, ph, idx, L.ishost, c);
        [m1, ~] = lqn_ph_moments(alpha, T);
        L.svcmeanByClass(k) = m1;
        law = phLaw(alpha, T);
        for s = 1:numel(L.qstations)
            model.nodes{L.qstations(s)}.setService(cls, law);
        end
        clientDelay.setService(cls, Exp.fitMean(max(GlobalConstants.FineTol, phDelayMean(self, lqn, ph, idx, c))));
    end

    for r = 1:size(L.openArrivals,1)
        k = L.openArrivals(r,1);
        tag = L.openArrivals(r,2);
        cls = model.classes{k};
        if tag > 0
            % entry arrival: the processor demand law of the entry is static
            L.svcmeanByClass(k) = ph.hostmean(tag);
            continue
        end
        cidx = -tag;
        eidx = lqn.callpair(cidx,2);
        L.svcmeanByClass(k) = ph.entrymean(eidx);
        law = phLaw(ph.entryalpha{eidx}, ph.entryT{eidx});
        for s = 1:numel(L.qstations)
            model.nodes{L.qstations(s)}.setService(cls, law);
        end
        aidx = lqn.callpair(cidx,1);
        rate = self.tput(aidx) * lqn.callproc_mean(cidx);
        if ~isfinite(rate) || rate <= GlobalConstants.FineTol
            rate = GlobalConstants.FineTol;
        end
        model.nodes{model.attribute.sourceIdx}.setArrival(cls, Exp.fitRate(rate));
    end

    ph.layer{idx} = L;
end

self.ph = ph;
end

% ------------------------------------------------------------------------
function [alpha, T] = phServiceLaw(self, lqn, ph, idx, ishost, c) %#ok<INUSL>
% Law of the demand caller C places on the server of layer IDX per invocation
if ishost
    % mixture over the entries of C, weighted by their share of its requests
    entries = lqn.entriesof{c};
    alphas = {}; Ts = {}; probs = [];
    for eidx = entries
        if isempty(ph.hostT{eidx}) || ph.share(eidx) <= 0
            continue
        end
        alphas{end+1} = ph.hostalpha{eidx}; %#ok<AGROW>
        Ts{end+1} = ph.hostT{eidx}; %#ok<AGROW>
        probs(end+1) = ph.share(eidx); %#ok<AGROW>
    end
    if isempty(alphas)
        alpha = 1; T = -GlobalConstants.Immediate;
        return
    end
    probs = probs / sum(probs);
    [alpha, T] = Workflow.composeMixture(alphas, Ts, probs);
    return
end

% task layer: the total demand is the sum, over the entries of the server, of a
% geometric compound of the entry law of mean equal to the number of calls
alpha = []; T = [];
for eidx = lqn.entriesof{idx}
    n = ph.ncalls(c, eidx);
    if n <= GlobalConstants.FineTol || isempty(ph.entryT{eidx})
        continue
    end
    [a2, T2] = Workflow.composeLoopGeometric(ph.entryalpha{eidx}, ph.entryT{eidx}, n);
    if isempty(alpha)
        alpha = a2; T = T2;
    else
        [alpha, T] = Workflow.composeSerial(alpha, T, a2, T2);
    end
end
if isempty(alpha)
    alpha = 1; T = -GlobalConstants.Immediate;
end
end

% ------------------------------------------------------------------------
function z = phDelayMean(self, lqn, ph, idx, c)
% Mean time a thread of caller C spends away from the stations of the model that
% holds server IDX, per invocation: idle, plus whatever of its cycle that model
% does not hold as a station of its own.
%
% This is ONE closure for both layerings. Under 'srvn.ph' the model holds a
% single server, so a host layer charges the whole call burst to the delay and a
% task layer charges the caller's processor plus every other callee. Under
% 'flat.ph' the model holds every server, and only the think times are left.
%
% Every term is SUMMED in rather than obtained by subtracting from a total. That
% subtraction cancels catastrophically once a call time is large: a caller whose
% only callee is this server has the two terms equal, and 7 + 1.4e47 - 1.4e47 is
% 0, not 7, because the think time falls below the ULP of the call time. The
% layer then sees a client delay of zero, saturates, reports a residence time
% that inflates the very call time that caused the cancellation, and the fixed
% point runs away -- lqn_sockshop reached RespT 1.4e47 this way.
z = self.thinkt(c) + lqn_ref_thinktime(lqn, c);
if ~isfinite(z) || z < 0
    z = 0;
end
z = z + ph.actthinkt(c);
% the caller's own processor residence, unless this model holds that processor
hidx = lqn.parent(c);
if ~phServedHere(self, idx, hidx)
    z = z + ph.procresid(c);
end
% and the time spent at every callee this model does not hold
for t = 1:lqn.ntasks
    tidx = lqn.tshift + t;
    if phServedHere(self, idx, tidx)
        continue
    end
    z = z + ph.calltime(c, tidx);
end
if ~isfinite(z) || z < 0
    z = GlobalConstants.FineTol;
end
end

% ------------------------------------------------------------------------
function tf = phServedHere(self, idx, elem)
% True when LQN element ELEM is a station of the same model that holds server
% IDX. Under 'srvn.ph' that is ELEM == IDX; under 'flat.ph' it is every server.
tf = false;
if isnan(elem) || elem < 1 || elem > numel(self.idxhash)
    return
end
tf = ~isnan(self.idxhash(elem)) && self.idxhash(elem) == self.idxhash(idx);
end

% ------------------------------------------------------------------------
function law = phLaw(alpha, T)
% Station law of a composed workflow. A geometric loop over a body of two or
% more phases closes a cycle in the phase graph, and a cyclic generator is a PH
% and not an APH: no layer solver declares PH, so such a law is reduced to the
% APH with the SAME first two moments. AMVA and NC read exactly those two, so
% the reduction is lossless for them and is a two-moment fit for the
% phase-aware layer solvers.
if Workflow.isAcyclicGenerator(T)
    law = APH(alpha, T);
    return
end
[m1, scv] = lqn_ph_moments(alpha, T);
law = APH.fitMeanAndSCV(m1, scv);
end
