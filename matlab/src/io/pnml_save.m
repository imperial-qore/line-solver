function pnml_save(model, filename)
% PNML_SAVE Write a LINE Petri net to a PNML (ISO/IEC 15909-2) file.
%
%   PNML_SAVE(MODEL, FILENAME) writes the Place/Transition net held by MODEL to
%   FILENAME in the PNML place/transition grammar
%   http://www.pnml.org/version-2009/grammar/ptnet, so that a LINE net can be
%   read by the tools built around that corpus (GreatSPN, TINA, the Model
%   Checking Contest harnesses).
%
%   THE P/T GRAMMAR IS UNCOLOURED, so what it can carry is narrower than what
%   LINE can express, and the difference is REFUSED rather than approximated:
%
%     * more than one job class -- a LINE net whose tokens carry a class is a
%       coloured net, and flattening the colours would change the model;
%     * an open class, a Source or a Sink -- the P/T grammar has no unbounded
%       token source;
%     * a queueing place, which is a station rather than a place;
%     * a firing-rate dependence, which no PNML element can carry;
%     * a distribution outside the scalar-parameter families listed in
%       PNML_DIST_TO_XML below, e.g. a phase-type or a Markovian arrival
%       process given by matrices.
%
%   TIMING RIDES IN A TOOLSPECIFIC BLOCK, which is where the grammar puts
%   anything it does not define. Each LINE MODE becomes one PNML transition, so
%   that the arcs of a mode are the arcs of a transition as the grammar
%   requires; the block records which LINE transition and mode the PNML
%   transition came from, so PNML_LOAD regroups the modes it split. A reader
%   that ignores the block still sees a correct untimed P/T net.
%
%   See also PNML_LOAD, LINEMODEL_SAVE.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

if ~isa(model, 'Network')
    line_error(mfilename, 'pnml_save expects a Network holding a Petri net.');
end

nodes = model.getNodes();
classes = model.getClasses();

if numel(classes) ~= 1
    line_error(mfilename, sprintf(['The PNML place/transition grammar is UNCOLOURED, so it cannot carry a net ' ...
        'with %d job classes: its tokens are indistinguishable. Export a single-class net, or use ' ...
        'linemodel_save for the full model.'], numel(classes)));
end
if isa(classes{1}, 'OpenClass')
    line_error(mfilename, ['The PNML place/transition grammar has no unbounded token source, so an open ' ...
        'class cannot be represented. Close the class, or use linemodel_save.']);
end

places = {};
transitions = {};
for i = 1:numel(nodes)
    nd = nodes{i};
    if isa(nd, 'Place')
        if ~isempty(nd.queueing) && nd.queueing
            line_error(mfilename, sprintf(['Place %s is a QUEUEING place, i.e. a station with an embedded ' ...
                'queue, which the place/transition grammar cannot represent.'], nd.name));
        end
        places{end+1} = nd; %#ok<AGROW>
    elseif isa(nd, 'Transition')
        transitions{end+1} = nd; %#ok<AGROW>
    else
        line_error(mfilename, sprintf(['Node %s is a %s. A PNML place/transition net holds only places and ' ...
            'transitions; a Source, a Sink or a queueing station has no counterpart in the grammar.'], ...
            nd.name, class(nd)));
    end
end
if isempty(places)
    line_error(mfilename, 'The model holds no Place, so there is no Petri net to write.');
end

placeIndex = zeros(1, numel(places));
for p = 1:numel(places)
    placeIndex(p) = model.getNodeIndex(places{p}.name);
end

marking = pnml_initial_marking(model, places, classes{1});

netName = model.getName();
if isempty(netName)
    netName = 'net';
end

% THE DOCUMENT IS BUILT WHOLE, THEN WRITTEN. Some refusals -- an unrepresentable
% distribution, a marking-dependent firing rate -- are only reachable while
% walking the modes, and streaming to the file would leave a truncated PNML
% behind on the way out. This is also what the Java, Python and C++ writers do.
sb = {};
sb{end+1} = '<?xml version="1.0" encoding="UTF-8"?>';
sb{end+1} = '<pnml xmlns="http://www.pnml.org/version-2009/grammar/pnml">';
sb{end+1} = sprintf('  <net id="%s" type="http://www.pnml.org/version-2009/grammar/ptnet">', pnml_escape(netName));
sb{end+1} = sprintf('    <name><text>%s</text></name>', pnml_escape(netName));
sb{end+1} = '    <page id="page0">';

for p = 1:numel(places)
    sb{end+1} = sprintf('      <place id="%s">', pnml_escape(places{p}.name)); %#ok<AGROW>
    sb{end+1} = sprintf('        <name><text>%s</text></name>', pnml_escape(places{p}.name)); %#ok<AGROW>
    sb{end+1} = sprintf('        <initialMarking><text>%d</text></initialMarking>', marking(p)); %#ok<AGROW>
    sb{end+1} = '      </place>'; %#ok<AGROW>
end

arcId = 0;
arcs = {};
for t = 1:numel(transitions)
    tr = transitions{t};
    nmodes = tr.getNumberOfModes();
    if nmodes == 0
        line_error(mfilename, sprintf('Transition %s declares no mode, so it has no firing behaviour to write.', tr.name));
    end
    for m = 1:nmodes
        if ~isempty(tr.firingRateDependence) && numel(tr.firingRateDependence) >= m && ~isempty(tr.firingRateDependence{m})
            line_error(mfilename, sprintf(['Transition %s mode %d declares a marking-dependent firing rate, which no ' ...
                'PNML element can carry.'], tr.name, m));
        end
        tid = pnml_mode_id(tr, m, nmodes);
        sb{end+1} = sprintf('      <transition id="%s">', pnml_escape(tid)); %#ok<AGROW>
        sb{end+1} = sprintf('        <name><text>%s</text></name>', pnml_escape(tid)); %#ok<AGROW>
        sb{end+1} = '        <toolspecific tool="LINE" version="3.0">'; %#ok<AGROW>
        timing = TimingStrategy.toText(tr.timingStrategies(m));
        sb{end+1} = sprintf('          <mode transition="%s" name="%s" timing="%s" servers="%s" priority="%s" weight="%s">', ...
            pnml_escape(tr.name), pnml_escape(pnml_mode_name(tr, m)), timing, ...
            pnml_num(tr.numberOfServers(m)), pnml_num(tr.firingPriorities(m)), pnml_num(tr.firingWeights(m))); %#ok<AGROW>
        if tr.timingStrategies(m) == TimingStrategy.TIMED
            sb = [sb, pnml_dist_to_xml(tr.distributions{m}, tr.name, m)]; %#ok<AGROW>
        end
        sb{end+1} = '          </mode>'; %#ok<AGROW>
        sb{end+1} = '        </toolspecific>'; %#ok<AGROW>
        sb{end+1} = '      </transition>'; %#ok<AGROW>

        for p = 1:numel(places)
            pidx = placeIndex(p);
            w = tr.enablingConditions{m}(pidx, 1);
            if w > 0
                arcId = arcId + 1;
                arcs = [arcs, pnml_arc(arcId, places{p}.name, tid, round(w), false)]; %#ok<AGROW>
            end
            inh = tr.inhibitingConditions{m}(pidx, 1);
            if isfinite(inh)
                arcId = arcId + 1;
                arcs = [arcs, pnml_arc(arcId, places{p}.name, tid, round(inh), true)]; %#ok<AGROW>
            end
            f = tr.firingOutcomes{m}(pidx, 1);
            if f > 0
                arcId = arcId + 1;
                arcs = [arcs, pnml_arc(arcId, tid, places{p}.name, round(f), false)]; %#ok<AGROW>
            end
        end
    end
end

sb = [sb, arcs];
sb{end+1} = '    </page>';
sb{end+1} = '  </net>';
sb{end+1} = '</pnml>';

fid = fopen(filename, 'w');
if fid == -1
    line_error(mfilename, sprintf('Cannot open %s for writing.', filename));
end
fprintf(fid, '%s\n', sb{:});
fclose(fid);
end

function lines = pnml_arc(arcId, source, target, weight, inhibitor)
% One arc. An inhibitor arc is not in the P/T grammar itself; <type
% value="inhibitor"/> is the extension GreatSPN, TINA and PIPE all read, so it is
% the one written here.
lines = {sprintf('      <arc id="a%d" source="%s" target="%s">', arcId, pnml_escape(source), pnml_escape(target))};
if inhibitor
    lines{end+1} = '        <type value="inhibitor"/>';
end
lines{end+1} = sprintf('        <inscription><text>%d</text></inscription>', weight);
lines{end+1} = '      </arc>';
end

function marking = pnml_initial_marking(model, places, jobclass)
% Token count of each place in the initial marking. The state set on the place
% is authoritative; a place with no state holds the class population when it is
% the reference station of the class, and nothing otherwise, which is the
% default LINE itself applies.
marking = zeros(1, numel(places));
refstat = jobclass.refstat;
for p = 1:numel(places)
    st = places{p}.getState();
    if ~isempty(st)
        marking(p) = round(st(1, 1));
    elseif ~isempty(refstat) && strcmp(refstat.name, places{p}.name)
        marking(p) = round(jobclass.population);
    end
end
if sum(marking) == 0
    line_warning(mfilename, 'The initial marking is empty: no transition of the written net is enabled.\n');
end
end

function id = pnml_mode_id(tr, m, nmodes)
% One PNML transition per LINE MODE, since a mode carries its own arc weights
% and the grammar gives a transition one inscription per arc. The name is left
% alone for the single-mode case, which is the shape of an imported P/T net.
if nmodes == 1
    id = tr.name;
else
    id = sprintf('%s.%s', tr.name, pnml_mode_name(tr, m));
end
end

function nm = pnml_mode_name(tr, m)
if numel(tr.modeNames) >= m && ~isempty(tr.modeNames{m})
    nm = char(tr.modeNames{m});
else
    nm = sprintf('Mode%d', m);
end
end

function s = pnml_num(v)
% Numbers are written in the shortest form that reads back exactly, so an
% integral count does not acquire a decimal point and Inf keeps the spelling
% PNML_LOAD expects.
if isinf(v)
    if v > 0
        s = 'Inf';
    else
        s = '-Inf';
    end
elseif v == round(v)
    s = sprintf('%d', round(v));
else
    s = sprintf('%.17g', v);
end
end

function lines = pnml_dist_to_xml(dist, trname, m)
% Timing of a timed mode, as a distribution name and its scalar parameters.
%
% The parameter NAMES are LINE's own (params{i}.paramName), so the block is
% self-describing and PNML_LOAD reconstructs the object by name rather than by
% position. A distribution whose parameters are not scalars -- a phase-type or a
% MAP given by matrices -- is REFUSED here: writing its mean rate instead would
% produce a file that reads back as a different model.
if isempty(dist)
    line_error(mfilename, sprintf('Transition %s mode %d is timed but carries no distribution.', trname, m));
end
if isa(dist, 'Immediate')
    lines = {'            <distribution name="Immediate"/>'};
    return
end
if isa(dist, 'Disabled')
    lines = {'            <distribution name="Disabled"/>'};
    return
end
params = dist.params;
if isempty(params)
    line_error(mfilename, sprintf(['Transition %s mode %d holds a %s, which declares no scalar parameter and ' ...
        'therefore cannot be written to PNML.'], trname, m, class(dist)));
end
lines = {sprintf('            <distribution name="%s">', pnml_escape(class(dist)))};
for i = 1:numel(params)
    v = params{i}.paramValue;
    if ~(isnumeric(v) && isscalar(v))
        line_error(mfilename, sprintf(['Transition %s mode %d holds a %s whose parameter "%s" is not a scalar. ' ...
            'The PNML timing block carries scalar parameters only; a matrix-parameterized law (PH, APH, MAP, ' ...
            'MMPP2, ME, RAP) has no PNML representation and writing its mean rate instead would read back as a ' ...
            'different model.'], trname, m, class(dist), params{i}.paramName));
    end
    lines{end+1} = sprintf('              <parameter name="%s" value="%s"/>', ...
        pnml_escape(params{i}.paramName), pnml_num(v)); %#ok<AGROW>
end
lines{end+1} = '            </distribution>';
end

function s = pnml_escape(s)
s = char(s);
s = strrep(s, '&', '&amp;');
s = strrep(s, '<', '&lt;');
s = strrep(s, '>', '&gt;');
s = strrep(s, '"', '&quot;');
end
