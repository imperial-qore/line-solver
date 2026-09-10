"""
PNML (ISO/IEC 15909-2) place/transition nets, read and written.

Native-Python twin of matlab/src/io/pnml_save.m and pnml_load.m and of
jar/src/main/java/jline/io/PnmlIO.java. The grammar written is
http://www.pnml.org/version-2009/grammar/ptnet, so that a LINE net can be read
by the tools built around that corpus (GreatSPN, TINA, the Model Checking
Contest harnesses) and a net from that corpus can be analysed here.

THE P/T GRAMMAR IS UNCOLOURED, so what it can carry is narrower than what LINE
can express, and the difference is REFUSED rather than approximated: more than
one job class, an open class or a Source/Sink, a queueing place, a firing-rate
dependence, and any distribution outside the scalar-parameter families listed in
_DIST_PARAMS.

TIMING RIDES IN A TOOLSPECIFIC BLOCK, which is where the grammar puts what it
does not define. Each LINE MODE becomes one PNML transition, so that the arcs of
a mode are the arcs of a transition as the grammar requires; the block records
which LINE transition and mode the PNML transition came from, so the reader
regroups the modes the writer split. A reader that ignores the block still sees
a correct untimed P/T net, and a P/T net with no such block is read with every
transition TIMED and EXPONENTIAL AT RATE 1, the convention of the stochastic
Petri net literature.

Key functions:
    save_pnml: Write the Petri net of a Network to a PNML file
    load_pnml: Read a PNML place/transition net into a Network
"""

import xml.etree.ElementTree as ET
from typing import Optional

import numpy as np

from ..distributions import (Det, Disabled, Erlang, Exp, Gamma, HyperExp,
                             Immediate, Lognormal, Pareto, Uniform, Weibull)
from ..lang.classes import ClosedClass, OpenClass
from ..lang.network import Network
from ..lang.nodes import Place, TimingStrategy, Transition


def _dist_params(dist):
    """
    Name and scalar parameters of a distribution, under the parameter names
    MATLAB stores, which are the names the file carries in every codebase.

    Returning None means the distribution has no PNML representation; the caller
    turns that into a refusal naming the law, because writing its mean rate
    instead would produce a file that reads back as a different model.
    """
    name = type(dist).__name__
    if isinstance(dist, Immediate):
        return 'Immediate', []
    if isinstance(dist, Disabled):
        return 'Disabled', []
    if isinstance(dist, Exp):
        return 'Exp', [('lambda', float(dist.get_rate()))]
    if isinstance(dist, Erlang):
        return 'Erlang', [('alpha', float(dist._phase_rate)), ('r', float(dist._phases))]
    if isinstance(dist, HyperExp):
        if len(dist._probs) != 2:
            return None, None
        return 'HyperExp', [('p', float(dist._probs[0])),
                            ('lambda1', float(dist._rates[0])),
                            ('lambda2', float(dist._rates[1]))]
    if isinstance(dist, Det):
        return 'Det', [('t', float(dist._value))]
    if isinstance(dist, Uniform):
        return 'Uniform', [('min', float(dist._min)), ('max', float(dist._max))]
    if isinstance(dist, Gamma):
        return 'Gamma', [('alpha', float(dist._shape)), ('beta', float(dist._scale))]
    if isinstance(dist, Pareto):
        return 'Pareto', [('alpha', float(dist._alpha)), ('k', float(dist._scale))]
    if isinstance(dist, Weibull):
        # LINE stores Weibull as (alpha=scale, r=shape) while the constructor
        # takes (shape, scale); naming both explicitly is what keeps a round trip
        # from transposing them.
        return 'Weibull', [('alpha', float(dist._scale)), ('r', float(dist._shape))]
    if isinstance(dist, Lognormal):
        return 'Lognormal', [('mu', float(dist._mu)), ('sigma', float(dist._sigma))]
    return None, name


def _dist_from(name: str, p: dict):
    """
    Rebuild the distribution of a timed mode. The constructors are named ONE BY
    ONE rather than applied positionally, for the Weibull reason above.
    """
    def need(key):
        if key not in p or p[key] is None or (isinstance(p[key], float) and np.isnan(p[key])):
            raise ValueError('Distribution %s is missing the parameter "%s".' % (name, key))
        return p[key]

    if name == 'Immediate':
        return Immediate()
    if name == 'Disabled':
        return Disabled()
    if name == 'Exp':
        return Exp(need('lambda'))
    if name == 'Det':
        return Det(need('t'))
    if name == 'Erlang':
        return Erlang(need('alpha'), int(round(need('r'))))
    if name == 'HyperExp':
        return HyperExp(need('p'), need('lambda1'), need('lambda2'))
    if name == 'Uniform':
        return Uniform(need('min'), need('max'))
    if name == 'Gamma':
        return Gamma(need('alpha'), need('beta'))
    if name == 'Pareto':
        return Pareto(need('alpha'), need('k'))
    if name == 'Weibull':
        return Weibull(need('r'), need('alpha'))
    if name == 'Lognormal':
        return Lognormal(need('mu'), need('sigma'))
    raise ValueError('The timing block names distribution "%s", which the PNML reader does not construct. '
                     'The families it reads are Exp, Det, Erlang, HyperExp, Uniform, Gamma, Pareto, Weibull, '
                     'Lognormal, Immediate and Disabled, which are the ones the writer writes.' % name)


def _num(v) -> str:
    """
    Shortest form that reads back exactly, so an integral count does not acquire
    a decimal point and an infinite server count keeps the spelling the reader
    expects.
    """
    v = float(v)
    if np.isinf(v):
        return 'Inf' if v > 0 else '-Inf'
    if v == round(v) and abs(v) < 9.007199254740992e15:
        return '%d' % int(round(v))
    return repr(v)


def _escape(s: str) -> str:
    return (str(s).replace('&', '&amp;').replace('<', '&lt;')
            .replace('>', '&gt;').replace('"', '&quot;'))


def _initial_marking(model, places, jobclass):
    """
    Token count of each place in the initial marking. The state set on the place
    is authoritative; a place with no state holds the class population when it is
    the reference station of the class, which is the default LINE itself applies.
    """
    marking = []
    refstat = getattr(jobclass, '_refstat', None)
    for pl in places:
        st = pl.get_state()
        if st is not None and np.size(st) > 0:
            marking.append(int(round(float(np.asarray(st).ravel()[0]))))
        elif refstat is not None and refstat.get_name() == pl.get_name():
            marking.append(int(round(float(jobclass.get_population()))))
        else:
            marking.append(0)
    return marking


def save_pnml(model, filename: str) -> None:
    """
    Write the Petri net held by a Network to a PNML place/transition file.

    Args:
        model: Network holding only places and transitions, with one closed class
        filename: output path

    Raises:
        ValueError: when the model holds something the grammar cannot carry
    """
    if not isinstance(model, Network):
        raise ValueError('save_pnml expects a Network holding a Petri net.')

    classes = model.get_classes()
    if len(classes) != 1:
        raise ValueError('The PNML place/transition grammar is UNCOLOURED, so it cannot carry a net with %d job '
                         'classes: its tokens are indistinguishable. Export a single-class net, or use save_model '
                         'for the full model.' % len(classes))
    if isinstance(classes[0], OpenClass):
        raise ValueError('The PNML place/transition grammar has no unbounded token source, so an open class cannot '
                         'be represented. Close the class, or use save_model.')

    places, transitions = [], []
    for nd in model.get_nodes():
        if isinstance(nd, Place):
            if getattr(nd, '_queueing', False):
                raise ValueError('Place %s is a QUEUEING place, i.e. a station with an embedded queue, which the '
                                 'place/transition grammar cannot represent.' % nd.get_name())
            places.append(nd)
        elif isinstance(nd, Transition):
            transitions.append(nd)
        else:
            raise ValueError('Node %s is a %s. A PNML place/transition net holds only places and transitions; a '
                             'Source, a Sink or a queueing station has no counterpart in the grammar.'
                             % (nd.get_name(), type(nd).__name__))
    if not places:
        raise ValueError('The model holds no Place, so there is no Petri net to write.')

    marking = _initial_marking(model, places, classes[0])
    net_name = model.get_name() or 'net'

    out = ['<?xml version="1.0" encoding="UTF-8"?>',
           '<pnml xmlns="http://www.pnml.org/version-2009/grammar/pnml">',
           '  <net id="%s" type="http://www.pnml.org/version-2009/grammar/ptnet">' % _escape(net_name),
           '    <name><text>%s</text></name>' % _escape(net_name),
           '    <page id="page0">']

    for i, pl in enumerate(places):
        out.append('      <place id="%s">' % _escape(pl.get_name()))
        out.append('        <name><text>%s</text></name>' % _escape(pl.get_name()))
        out.append('        <initialMarking><text>%d</text></initialMarking>' % marking[i])
        out.append('      </place>')

    arcs = []
    arc_id = 0
    for tr in transitions:
        nmodes = tr.get_number_of_modes()
        if nmodes == 0:
            raise ValueError('Transition %s declares no mode, so it has no firing behaviour to write.' % tr.get_name())
        for m in range(nmodes):
            if tr._firing_rate_dependence[m] is not None:
                raise ValueError('Transition %s mode %d declares a marking-dependent firing rate, which no PNML '
                                 'element can carry.' % (tr.get_name(), m + 1))
            mode_name = tr._mode_names[m] or ('Mode%d' % (m + 1))
            tid = tr.get_name() if nmodes == 1 else '%s.%s' % (tr.get_name(), mode_name)
            timing = tr._timing_strategies[m]
            out.append('      <transition id="%s">' % _escape(tid))
            out.append('        <name><text>%s</text></name>' % _escape(tid))
            out.append('        <toolspecific tool="LINE" version="3.0">')
            out.append('          <mode transition="%s" name="%s" timing="%s" servers="%s" priority="%s" weight="%s">'
                       % (_escape(tr.get_name()), _escape(mode_name),
                          'immediate' if timing == TimingStrategy.IMMEDIATE else 'timed',
                          _num(tr._number_of_servers[m]), _num(tr._firing_priorities[m]),
                          _num(tr._firing_weights[m])))
            if timing != TimingStrategy.IMMEDIATE:
                dist = tr._distributions[m]
                if dist is None:
                    raise ValueError('Transition %s mode %d is timed but carries no distribution.'
                                     % (tr.get_name(), m + 1))
                name, params = _dist_params(dist)
                if name is None:
                    raise ValueError('Transition %s mode %d holds a %s, whose parameters are not scalars. The PNML '
                                     'timing block carries scalar parameters only; a matrix-parameterized law (PH, '
                                     'APH, MAP, MMPP2, ME, RAP) has no PNML representation and writing its mean rate '
                                     'instead would read back as a different model.'
                                     % (tr.get_name(), m + 1, params))
                if not params:
                    out.append('            <distribution name="%s"/>' % _escape(name))
                else:
                    out.append('            <distribution name="%s">' % _escape(name))
                    for key, value in params:
                        out.append('              <parameter name="%s" value="%s"/>' % (_escape(key), _num(value)))
                    out.append('            </distribution>')
            out.append('          </mode>')
            out.append('        </toolspecific>')
            out.append('      </transition>')

            for p, pl in enumerate(places):
                # get_node_index is 1-BASED while the mode matrices are 0-based
                # numpy arrays, so the row is one less than the node index.
                pidx = model.get_node_index(pl) - 1
                w = tr._enabling_conditions[m][pidx, 0]
                if w > 0:
                    arc_id += 1
                    arcs.append(_arc(arc_id, pl.get_name(), tid, int(round(w)), False))
                inh = tr._inhibiting_conditions[m][pidx, 0]
                if np.isfinite(inh):
                    arc_id += 1
                    arcs.append(_arc(arc_id, pl.get_name(), tid, int(round(inh)), True))
                f = tr._firing_outcomes[m][pidx, 0]
                if f > 0:
                    arc_id += 1
                    arcs.append(_arc(arc_id, tid, pl.get_name(), int(round(f)), False))

    out.extend(arcs)
    out.append('    </page>')
    out.append('  </net>')
    out.append('</pnml>')

    with open(filename, 'w', encoding='utf-8') as fh:
        fh.write('\n'.join(out) + '\n')


def _arc(arc_id: int, source: str, target: str, weight: int, inhibitor: bool):
    lines = ['      <arc id="a%d" source="%s" target="%s">' % (arc_id, _escape(source), _escape(target))]
    if inhibitor:
        # An inhibitor arc is not in the P/T grammar itself; <type value="inhibitor"/>
        # is the extension GreatSPN, TINA and PIPE all read, so it is the one written.
        lines.append('        <type value="inhibitor"/>')
    lines.append('        <inscription><text>%d</text></inscription>' % weight)
    lines.append('      </arc>')
    return '\n'.join(lines)


def _local(tag: str) -> str:
    return tag.split('}')[-1].split(':')[-1]


def _children(elem, local_name: str):
    return [c for c in list(elem) if _local(c.tag) == local_name]


def _descendants(elem, local_name: str):
    out = []
    for c in list(elem):
        if _local(c.tag) == local_name:
            out.append(c)
        else:
            out.extend(_descendants(c, local_name))
    return out


def _text_of(elem) -> str:
    texts = _children(elem, 'text')
    if texts:
        return (texts[0].text or '').strip()
    return ''.join(elem.itertext()).strip()


def _text_number(elem, label: str, default: float) -> float:
    labels = _children(elem, label)
    if not labels:
        return default
    txt = _text_of(labels[0])
    if not txt:
        return default
    # Some tools write the marking as "3" and some as "1`3" (a coloured multiset
    # of one colour); the plain integer is the one the grammar defines.
    if '`' in txt:
        txt = txt.rsplit('`', 1)[1]
    try:
        return float(txt.strip())
    except ValueError:
        raise ValueError('Label <%s> holds "%s", which is not a number.' % (label, txt))


def _attr_number(elem, key: str, default: float) -> float:
    txt = (elem.get(key) or '').strip()
    if not txt:
        return default
    if txt.lower() == 'inf':
        return float('inf')
    if txt.lower() == '-inf':
        return float('-inf')
    try:
        return float(txt)
    except ValueError:
        raise ValueError('Attribute %s holds "%s", which is not a number.' % (key, txt))


def _element_id(elem) -> str:
    eid = elem.get('id') or ''
    if not eid:
        names = _children(elem, 'name')
        if names:
            eid = _text_of(names[0])
    if not eid:
        raise ValueError('A place or transition carries neither an id nor a name.')
    return eid


def _is_inhibitor(arc_elem) -> bool:
    for t in _children(arc_elem, 'type'):
        if (t.get('value') or '').lower() == 'inhibitor':
            return True
    return (arc_elem.get('type') or '').lower() == 'inhibitor'


def _line_toolspecific(trans_elem):
    for block in _children(trans_elem, 'toolspecific'):
        if (block.get('tool') or '').upper() != 'LINE':
            continue
        modes = _children(block, 'mode')
        if modes:
            return modes[0]
    return None


def load_pnml(filename: str, net_id: Optional[str] = None) -> Network:
    """
    Read one net of a PNML place/transition document into a LINE Network.

    Args:
        filename: input path
        net_id: id of the net to read, or None for the first

    Returns:
        The equivalent Network: one Place per PNML place, one Transition per
        group of PNML transitions belonging to the same LINE transition, and a
        single ClosedClass holding the tokens of the initial marking.
    """
    root = ET.parse(filename).getroot()
    nets = _children(root, 'net')
    if not nets:
        raise ValueError('%s holds no <net> element.' % filename)
    net = None
    for candidate in nets:
        if not net_id or candidate.get('id') == net_id:
            net = candidate
            break
    if net is None:
        raise ValueError('%s holds no net with id "%s".' % (filename, net_id))
    nettype = net.get('type') or ''
    if nettype and 'ptnet' not in nettype:
        raise ValueError('Net "%s" declares type %s. Only the place/transition grammar '
                         '(http://www.pnml.org/version-2009/grammar/ptnet) is read: a coloured or a symmetric net '
                         'carries token colours that a single-class LINE net cannot hold.' % (net.get('id'), nettype))

    # Places, transitions and arcs may sit directly under <net> or under any
    # <page>; the grammar allows both and tools differ, so the whole subtree is
    # searched rather than one level of it.
    place_elems = _descendants(net, 'place')
    trans_elems = _descendants(net, 'transition')
    arc_elems = _descendants(net, 'arc')
    if not place_elems:
        raise ValueError('Net "%s" holds no place.' % net.get('id'))

    place_names = [_element_id(pe) for pe in place_elems]
    place_marking = [int(_text_number(pe, 'initialMarking', 0.0)) for pe in place_elems]
    if sum(place_marking) == 0:
        raise ValueError('The initial marking of this net is empty. A LINE closed class needs tokens to hold, and a '
                         'net with none has no reachable behaviour to analyse.')

    trans_ids, mode_owner, mode_name = [], [], []
    mode_timing, mode_servers, mode_priority, mode_weight, mode_dist = [], [], [], [], []
    for te in trans_elems:
        tid = _element_id(te)
        trans_ids.append(tid)
        owner, name = tid, 'Mode1'
        timing, servers, priority, weight = TimingStrategy.TIMED, 1, 1.0, 1.0
        dist = Exp(1)
        spec = _line_toolspecific(te)
        if spec is not None:
            owner = spec.get('transition') or owner
            name = spec.get('name') or name
            if (spec.get('timing') or '').lower() == 'immediate':
                timing = TimingStrategy.IMMEDIATE
            sv = _attr_number(spec, 'servers', 1.0)
            servers = int(round(sv)) if np.isfinite(sv) else int(2 ** 31 - 1)
            priority = _attr_number(spec, 'priority', 1.0)
            weight = _attr_number(spec, 'weight', 1.0)
            if timing != TimingStrategy.IMMEDIATE:
                dists = _children(spec, 'distribution')
                if dists:
                    params = {}
                    for pe in _children(dists[0], 'parameter'):
                        params[pe.get('name')] = _attr_number(pe, 'value', float('nan'))
                    dist = _dist_from(dists[0].get('name') or 'Exp', params)
        mode_owner.append(owner)
        mode_name.append(name)
        mode_timing.append(timing)
        mode_servers.append(servers)
        mode_priority.append(priority)
        mode_weight.append(weight)
        mode_dist.append(dist)

    owner_names, owner_of = [], []
    for owner in mode_owner:
        if owner not in owner_names:
            owner_names.append(owner)
        owner_of.append(owner_names.index(owner))

    name = net_id or net.get('id') or 'pnml'
    model = Network(name)
    place_obj = [Place(model, nm) for nm in place_names]
    trans_obj = [Transition(model, nm) for nm in owner_names]

    # The reference station is the first place holding tokens, so the class
    # starts where the marking says it does.
    ref_idx = next(i for i, mk in enumerate(place_marking) if mk > 0)
    jobclass = ClosedClass(model, 'Class1', int(sum(place_marking)), place_obj[ref_idx], 0)

    mode_obj = []
    for i in range(len(trans_ids)):
        tr = trans_obj[owner_of[i]]
        mode = tr.add_mode(mode_name[i])
        mode_obj.append(mode)
        tr.set_timing_strategy(mode, mode_timing[i])
        tr.set_number_of_servers(mode, mode_servers[i])
        tr.set_firing_priorities(mode, mode_priority[i])
        tr.set_firing_weights(mode, mode_weight[i])
        if mode_timing[i] != TimingStrategy.IMMEDIATE:
            tr.set_distribution(mode, mode_dist[i])

    routing = model.init_routing_matrix()
    for ae in arc_elems:
        src, tgt = ae.get('source'), ae.get('target')
        w = _text_number(ae, 'inscription', 1.0)
        if w <= 0:
            raise ValueError('Arc %s -> %s carries a non-positive inscription %g.' % (src, tgt, w))
        if src in place_names and tgt in trans_ids:
            ip, it = place_names.index(src), trans_ids.index(tgt)
            tr = trans_obj[owner_of[it]]
            if _is_inhibitor(ae):
                tr.set_inhibiting_conditions(mode_obj[it], jobclass, place_obj[ip], int(round(w)))
            else:
                tr.set_enabling_conditions(mode_obj[it], jobclass, place_obj[ip], int(round(w)))
            routing.set(jobclass, jobclass, place_obj[ip], tr, 1.0)
            continue
        if src in trans_ids and tgt in place_names:
            it, ip = trans_ids.index(src), place_names.index(tgt)
            tr = trans_obj[owner_of[it]]
            tr.set_firing_outcome(mode_obj[it], jobclass, place_obj[ip], int(round(w)))
            routing.set(jobclass, jobclass, tr, place_obj[ip], 1.0)
            continue
        raise ValueError('Arc %s -> %s connects two places or two transitions, which the place/transition grammar '
                         'does not allow.' % (src, tgt))
    model.link(routing)

    for i, pl in enumerate(place_obj):
        pl.set_state(place_marking[i])
    return model
