function model = pnml_load(filename, netName)
% PNML_LOAD Read a PNML (ISO/IEC 15909-2) place/transition net into LINE.
%
%   MODEL = PNML_LOAD(FILENAME) reads the first net of a PNML document in the
%   place/transition grammar and returns the equivalent LINE Network: one Place
%   per PNML place, one Transition per PNML transition, and a single ClosedClass
%   holding the tokens of the initial marking.
%
%   MODEL = PNML_LOAD(FILENAME, NETNAME) reads the net with the given id, for a
%   document holding more than one.
%
%   THE GRAMMAR IS UNTIMED. A PNML place/transition net says nothing about how
%   long a transition takes to fire, so a file written by another tool is read
%   with every transition TIMED and EXPONENTIAL AT RATE 1, which is the
%   convention of the stochastic Petri net literature and of GreatSPN's own
%   default. A file written by PNML_SAVE carries LINE's timing in its
%   toolspecific block and reads back with the distributions, servers,
%   priorities and weights it was written with, and with the modes regrouped
%   into the transitions they were split from.
%
%   Inhibitor arcs are read from the <type value="inhibitor"/> extension that
%   GreatSPN, TINA and PIPE all write, since the grammar itself has no
%   inhibitor arc.
%
%   See also PNML_SAVE, LINEMODEL_LOAD.
%
% Copyright (c) 2012-2026, Imperial College London
% All rights reserved.

import javax.xml.parsers.DocumentBuilderFactory; %#ok<JAPIEXT629>

if nargin < 2
    netName = '';
end

fid = fopen(filename, 'r');
if fid == -1
    line_error(mfilename, sprintf('File %s cannot be found.', filename));
end
fclose(fid);

dbFactory = DocumentBuilderFactory.newInstance();
dBuilder = dbFactory.newDocumentBuilder();
doc = dBuilder.parse(filename);

nets = pnml_children(doc.getDocumentElement(), 'net');
if isempty(nets)
    line_error(mfilename, sprintf('%s holds no <net> element.', filename));
end
net = [];
for i = 1:numel(nets)
    if isempty(netName) || strcmp(char(nets{i}.getAttribute('id')), netName)
        net = nets{i};
        break
    end
end
if isempty(net)
    line_error(mfilename, sprintf('%s holds no net with id "%s".', filename, netName));
end
nettype = char(net.getAttribute('type'));
if ~isempty(nettype) && isempty(strfind(nettype, 'ptnet')) %#ok<STREMP>
    line_error(mfilename, sprintf(['Net "%s" declares type %s. Only the place/transition grammar ' ...
        '(http://www.pnml.org/version-2009/grammar/ptnet) is read: a coloured or a symmetric net carries ' ...
        'token colours that a single-class LINE net cannot hold.'], char(net.getAttribute('id')), nettype));
end

% Places, transitions and arcs may sit directly under <net> or under any
% <page>; the grammar allows both and tools differ, so the whole subtree is
% searched rather than one level of it.
placeElems = pnml_descendants(net, 'place');
transElems = pnml_descendants(net, 'transition');
arcElems = pnml_descendants(net, 'arc');

if isempty(placeElems)
    line_error(mfilename, sprintf('Net "%s" holds no place.', char(net.getAttribute('id'))));
end

placeNames = cell(1, numel(placeElems));
placeMarking = zeros(1, numel(placeElems));
for i = 1:numel(placeElems)
    placeNames{i} = pnml_id(placeElems{i});
    placeMarking(i) = pnml_text_number(placeElems{i}, 'initialMarking', 0);
end
if sum(placeMarking) == 0
    line_error(mfilename, ['The initial marking of this net is empty. A LINE closed class needs tokens to ' ...
        'hold, and a net with none has no reachable behaviour to analyse.']);
end

% Each PNML transition is one LINE MODE. The toolspecific block, when present,
% says which LINE transition the mode belongs to, so the modes PNML_SAVE split
% into separate transitions are regrouped here.
nt = numel(transElems);
modeOwner = cell(1, nt);
modeName = cell(1, nt);
modeTiming = zeros(1, nt);
modeServers = ones(1, nt);
modePriority = ones(1, nt);
modeWeight = ones(1, nt);
modeDist = cell(1, nt);
transId = cell(1, nt);
for i = 1:nt
    transId{i} = pnml_id(transElems{i});
    modeOwner{i} = transId{i};
    modeName{i} = 'Mode1';
    modeTiming(i) = TimingStrategy.TIMED;
    modeDist{i} = Exp(1);
    spec = pnml_line_toolspecific(transElems{i});
    if ~isempty(spec)
        if ~isempty(spec.transition)
            modeOwner{i} = spec.transition;
        end
        if ~isempty(spec.name)
            modeName{i} = spec.name;
        end
        modeTiming(i) = spec.timing;
        modeServers(i) = spec.servers;
        modePriority(i) = spec.priority;
        modeWeight(i) = spec.weight;
        if modeTiming(i) == TimingStrategy.TIMED
            modeDist{i} = spec.distribution;
        end
    end
end

[ownerNames, ~, ownerOf] = unique(modeOwner, 'stable');

if isempty(netName)
    netName = char(net.getAttribute('id'));
end
if isempty(netName)
    netName = 'pnml';
end

model = Network(netName);

placeObj = cell(1, numel(placeNames));
for i = 1:numel(placeNames)
    placeObj{i} = Place(model, placeNames{i});
end
transObj = cell(1, numel(ownerNames));
for i = 1:numel(ownerNames)
    transObj{i} = Transition(model, ownerNames{i});
end

% The reference station is the first place holding tokens, so the class starts
% where the marking says it does.
refIdx = find(placeMarking > 0, 1, 'first');
jobclass = ClosedClass(model, 'Class1', sum(placeMarking), placeObj{refIdx}, 0);

% Modes, in the order the transitions appear in the file, so that a
% write/read/write cycle reproduces the file it started from.
modeIndex = zeros(1, nt);
for i = 1:nt
    tr = transObj{ownerOf(i)};
    tr.addMode(modeName{i});
    % addMode returns a Mode object, whose subsindex is 0-based; the setters
    % below index plain arrays, so the 1-based count is what is kept.
    modeIndex(i) = tr.getNumberOfModes();
    tr.setTimingStrategy(modeIndex(i), modeTiming(i));
    tr.setNumberOfServers(modeIndex(i), modeServers(i));
    tr.setFiringPriorities(modeIndex(i), modePriority(i));
    tr.setFiringWeights(modeIndex(i), modeWeight(i));
    if modeTiming(i) == TimingStrategy.TIMED
        tr.setDistribution(modeIndex(i), modeDist{i});
    end
end

placeByName = pnml_name_lookup(placeNames);
transByName = pnml_name_lookup(transId);

R = model.initRoutingMatrix();
for a = 1:numel(arcElems)
    src = char(arcElems{a}.getAttribute('source'));
    tgt = char(arcElems{a}.getAttribute('target'));
    w = pnml_text_number(arcElems{a}, 'inscription', 1);
    if w <= 0
        line_error(mfilename, sprintf('Arc %s -> %s carries a non-positive inscription %g.', src, tgt, w));
    end
    ip = placeByName(src);
    it = transByName(tgt);
    if ip > 0 && it > 0
        tr = transObj{ownerOf(it)};
        if pnml_is_inhibitor(arcElems{a})
            tr.setInhibitingConditions(modeIndex(it), jobclass, placeObj{ip}, w);
        else
            tr.setEnablingConditions(modeIndex(it), jobclass, placeObj{ip}, w);
        end
        R{1,1}(placeObj{ip}, tr) = 1;
        continue
    end
    it = transByName(src);
    ip = placeByName(tgt);
    if it > 0 && ip > 0
        tr = transObj{ownerOf(it)};
        tr.setFiringOutcome(modeIndex(it), jobclass, placeObj{ip}, w);
        R{1,1}(tr, placeObj{ip}) = 1;
        continue
    end
    line_error(mfilename, sprintf(['Arc %s -> %s connects two places or two transitions, which the ' ...
        'place/transition grammar does not allow.'], src, tgt));
end

model.link(R);

for i = 1:numel(placeObj)
    placeObj{i}.setState(placeMarking(i));
end
end

function idx = pnml_name_lookup(names)
% Name -> index lookup returning 0 for an absent key, as an anonymous function
% over a cell array of names. containers.Map is not used anywhere in LINE.
idx = @(key) pnml_index_of(names, key);
end

function k = pnml_index_of(names, key)
k = 0;
for i = 1:numel(names)
    if strcmp(names{i}, key)
        k = i;
        return
    end
end
end

function out = pnml_children(elem, localName)
% Direct children with the given local name, ignoring any namespace prefix.
out = {};
kids = elem.getChildNodes();
for i = 0:kids.getLength()-1
    nd = kids.item(i);
    if nd.getNodeType() == 1 && strcmp(pnml_local(char(nd.getNodeName())), localName)
        out{end+1} = nd; %#ok<AGROW>
    end
end
end

function out = pnml_descendants(elem, localName)
% Descendants with the given local name, in document order.
out = {};
kids = elem.getChildNodes();
for i = 0:kids.getLength()-1
    nd = kids.item(i);
    if nd.getNodeType() ~= 1
        continue
    end
    if strcmp(pnml_local(char(nd.getNodeName())), localName)
        out{end+1} = nd; %#ok<AGROW>
    else
        sub = pnml_descendants(nd, localName);
        for j = 1:numel(sub)
            out{end+1} = sub{j}; %#ok<AGROW>
        end
    end
end
end

function nm = pnml_local(nm)
k = strfind(nm, ':');
if ~isempty(k)
    nm = nm(k(end)+1:end);
end
end

function id = pnml_id(elem)
id = char(elem.getAttribute('id'));
if isempty(id)
    nameElems = pnml_children(elem, 'name');
    if ~isempty(nameElems)
        id = pnml_text_of(nameElems{1});
    end
end
if isempty(id)
    line_error(mfilename, 'A place or transition carries neither an id nor a name.');
end
end

function txt = pnml_text_of(elem)
% The <text> child of a label element, which is where the grammar puts the
% value; an element carrying bare character data is also accepted.
textElems = pnml_children(elem, 'text');
if ~isempty(textElems)
    txt = strtrim(char(textElems{1}.getTextContent()));
else
    txt = strtrim(char(elem.getTextContent()));
end
end

function v = pnml_text_number(elem, labelName, defaultValue)
v = defaultValue;
labels = pnml_children(elem, labelName);
if isempty(labels)
    return
end
txt = pnml_text_of(labels{1});
if isempty(txt)
    return
end
% Some tools write the marking as "3" and some as "1`3" (a coloured multiset of
% one colour); the plain integer is the one the grammar defines.
tick = strfind(txt, '`');
if ~isempty(tick)
    txt = txt(tick(end)+1:end);
end
num = str2double(txt);
if isnan(num)
    line_error(mfilename, sprintf('Label <%s> holds "%s", which is not a number.', labelName, txt));
end
v = num;
end

function tf = pnml_is_inhibitor(arcElem)
tf = false;
types = pnml_children(arcElem, 'type');
for i = 1:numel(types)
    if strcmpi(char(types{i}.getAttribute('value')), 'inhibitor')
        tf = true;
        return
    end
end
if strcmpi(char(arcElem.getAttribute('type')), 'inhibitor')
    tf = true;
end
end

function spec = pnml_line_toolspecific(transElem)
% The <toolspecific tool="LINE"> block of a transition, or empty when the file
% comes from another tool.
spec = [];
blocks = pnml_children(transElem, 'toolspecific');
for i = 1:numel(blocks)
    if ~strcmpi(char(blocks{i}.getAttribute('tool')), 'LINE')
        continue
    end
    modes = pnml_children(blocks{i}, 'mode');
    if isempty(modes)
        continue
    end
    md = modes{1};
    spec.transition = char(md.getAttribute('transition'));
    spec.name = char(md.getAttribute('name'));
    timing = char(md.getAttribute('timing'));
    if strcmpi(timing, 'immediate')
        spec.timing = TimingStrategy.IMMEDIATE;
    else
        spec.timing = TimingStrategy.TIMED;
    end
    spec.servers = pnml_attr_number(md, 'servers', 1);
    spec.priority = pnml_attr_number(md, 'priority', 1);
    spec.weight = pnml_attr_number(md, 'weight', 1);
    spec.distribution = pnml_dist_from_xml(md);
    return
end
end

function v = pnml_attr_number(elem, key, defaultValue)
txt = strtrim(char(elem.getAttribute(key)));
if isempty(txt)
    v = defaultValue;
    return
end
if strcmpi(txt, 'Inf')
    v = Inf;
    return
end
if strcmpi(txt, '-Inf')
    v = -Inf;
    return
end
v = str2double(txt);
if isnan(v)
    line_error(mfilename, sprintf('Attribute %s holds "%s", which is not a number.', key, txt));
end
end

function dist = pnml_dist_from_xml(modeElem)
% Rebuild the distribution of a timed mode from the toolspecific block.
%
% The constructors are named ONE BY ONE rather than applied positionally: LINE
% stores Weibull as (alpha=scale, r=shape) while its constructor takes
% (shape, scale), so a positional rebuild would silently transpose the two.
dist = Exp(1);
dists = pnml_children(modeElem, 'distribution');
if isempty(dists)
    return
end
d = dists{1};
name = char(d.getAttribute('name'));
p = struct();
params = pnml_children(d, 'parameter');
for i = 1:numel(params)
    key = char(params{i}.getAttribute('name'));
    p.(matlab.lang.makeValidName(key)) = pnml_attr_number(params{i}, 'value', NaN);
end

switch name
    case 'Immediate'
        dist = Immediate.getInstance();
    case 'Disabled'
        dist = Disabled.getInstance();
    case 'Exp'
        dist = Exp(p.lambda);
    case 'Det'
        dist = Det(p.t);
    case 'Erlang'
        dist = Erlang(p.alpha, p.r);
    case 'HyperExp'
        dist = HyperExp(p.p, p.lambda1, p.lambda2);
    case 'Uniform'
        dist = Uniform(p.min, p.max);
    case 'Gamma'
        dist = Gamma(p.alpha, p.beta);
    case 'Pareto'
        dist = Pareto(p.alpha, p.k);
    case 'Weibull'
        dist = Weibull(p.r, p.alpha);
    case 'Lognormal'
        dist = Lognormal(p.mu, p.sigma);
    otherwise
        line_error(mfilename, sprintf(['The timing block names distribution "%s", which PNML_LOAD does not ' ...
            'construct. The families it reads are Exp, Det, Erlang, HyperExp, Uniform, Gamma, Pareto, Weibull, ' ...
            'Lognormal, Immediate and Disabled, which are the ones PNML_SAVE writes.'], name));
end
end
