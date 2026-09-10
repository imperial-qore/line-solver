"""Two structural refusals that belong to getStruct, not to the parser.

Both documents are valid against `lqn-core.xsd` and both are accepted by lqns,
so nothing in the reader stops them; what they describe is a model LINE's
`LayeredNetworkStruct` cannot represent, and the refusal is by name.

  1. An entry whose `<entry-phase-activities>` is EMPTY reaches no activity, so
     it has no service and no reply. Python and the JAR used to build the
     struct anyway and report a row of NaN for it; MATLAB and C++ refused.
  2. An activity that REPLIES and then continues into a phase-1 successor ends
     phase 1 of its entry, so the tail is a second phase the struct has no slot
     for. MATLAB, the JAR and python refused; C++ had no guard.

The wording is `getStruct.m`'s, verbatim, in all four codebases. The twins are
cpp/tests/test_lqn_struct_guards.cpp, jar LqnStructGuardsTest.java and MATLAB
test_lqn_struct_guards.m. `user-models/blsr.lqnx` in the LQNS corpus carries
BOTH defects, and refuses on the first, as MATLAB does.
"""
import os
import tempfile

import pytest

from line_solver import LayeredNetwork

HEAD = ('<?xml version="1.0"?>\n'
        '<lqn-model name="g" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n')
TAIL = '</lqn-model>\n'

CLIENT = ('  <processor name="PC" scheduling="inf">\n'
          '    <task name="TC" scheduling="ref" multiplicity="2">\n'
          '      <entry name="EC" type="PH1PH2">\n'
          '        <entry-phase-activities>\n'
          '          <activity name="EC_ph1" phase="1" host-demand-mean="1">\n'
          '            <synch-call dest="ES" calls-mean="1"/>\n'
          '          </activity>\n'
          '        </entry-phase-activities>\n'
          '      </entry>\n'
          '    </task>\n'
          '  </processor>\n')

SERVER = ('    <task name="TS" scheduling="fcfs">\n'
          '      <entry name="ES" type="PH1PH2">\n'
          '        <entry-phase-activities>\n'
          '          <activity name="ES_ph1" phase="1" host-demand-mean="2"/>\n'
          '        </entry-phase-activities>\n'
          '      </entry>\n'
          '    </task>\n')

# an entry that names no activity at all
EMPTY_ENTRY = ('    <task name="TN" scheduling="inf">\n'
               '      <entry name="EN" type="PH1PH2">\n'
               '        <entry-phase-activities>\n'
               '        </entry-phase-activities>\n'
               '      </entry>\n'
               '    </task>\n')

# a1 replies, a2 follows it and carries no phase attribute, so it is phase 1
REPLY_THEN_CONTINUE = (
    '    <task name="TR" scheduling="fcfs">\n'
    '      <entry name="ER" type="NONE"/>\n'
    '      <task-activities>\n'
    '        <activity name="a1" bound-to-entry="ER" host-demand-mean="1"/>\n'
    '        <activity name="a2" host-demand-mean="1"/>\n'
    '        <precedence>\n'
    '          <pre><activity name="a1"/></pre>\n'
    '          <post><activity name="a2"/></post>\n'
    '        </precedence>\n'
    '        <reply-entry name="ER"><reply-activity name="a1"/></reply-entry>\n'
    '      </task-activities>\n'
    '    </task>\n')


def build(path, extra_task, server_call='ES'):
    body = HEAD + CLIENT.replace('dest="ES"', 'dest="%s"' % server_call)
    body += '  <processor name="PS" scheduling="fcfs">\n' + SERVER + extra_task
    body += '  </processor>\n' + TAIL
    with open(path, 'w') as fd:
        fd.write(body)
    return path


def refusal(path):
    """The message getStruct raises, or the empty string."""
    try:
        LayeredNetwork.parseXML(path).getStruct()
    except Exception as e:  # ValueError here, RuntimeException in the JAR
        return str(e)
    return ''


def test_entry_with_no_bound_activity_is_refused():
    d = tempfile.mkdtemp(prefix='lqn_guards_')
    path = build(os.path.join(d, 'empty_entry.lqnx'), EMPTY_ENTRY)
    assert refusal(path) == 'An entry does not have any boundTo activity.'


def test_reply_then_continue_is_refused():
    d = tempfile.mkdtemp(prefix='lqn_guards_')
    path = build(os.path.join(d, 'reply_tail.lqnx'), REPLY_THEN_CONTINUE, server_call='ER')
    msg = refusal(path)
    assert 'Unsupported replyTo in non-terminal activity' in msg


def test_a_sound_model_still_builds():
    d = tempfile.mkdtemp(prefix='lqn_guards_')
    path = build(os.path.join(d, 'ok.lqnx'), '')
    lqn = LayeredNetwork.parseXML(path).getStruct()
    assert lqn.nentries == 2
    assert lqn.ntasks == 2
