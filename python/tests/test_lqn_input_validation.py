"""
The .lqnx reader must refuse a structurally inconsistent document.

A defective input used to be accepted in silence and to surface much later, as
a NaN metric, an aliased element or a model with no customers at all. The reader
now names the defect at its source. The same fourteen documents below are
refused with the SAME message by the MATLAB, JAR and C++ readers, so this file
also pins the wording that keeps the four codebases interchangeable.
"""

import pytest

from line_solver import LayeredNetwork

HEAD = ('<?xml version="1.0"?>\n'
        '<lqn-model name="t" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n')
TAIL = '</lqn-model>\n'


def entry(name, calls=(), fwd=(), arrival=None):
    s = '      <entry name="%s" type="PH1PH2"' % name
    if arrival is not None:
        s += ' open-arrival-rate="%s"' % arrival
    s += '>\n'
    for dest, prob in fwd:
        s += '        <forwarding dest="%s" prob="%s"/>\n' % (dest, prob)
    s += '        <entry-phase-activities>\n'
    s += '          <activity name="%s_ph1" phase="1" host-demand-mean="1">\n' % name
    for dest in calls:
        s += '            <synch-call dest="%s" calls-mean="1"/>\n' % dest
    s += '          </activity>\n        </entry-phase-activities>\n      </entry>\n'
    return s


def task(name, sched, entries, extra=''):
    return '    <task name="%s" scheduling="%s">\n%s%s    </task>\n' % (
        name, sched, ''.join(entries), extra)


def proc(name, tasks, sched='fcfs'):
    return '  <processor name="%s" scheduling="%s">\n%s  </processor>\n' % (
        name, sched, ''.join(tasks))


REF = task('t0', 'ref', [entry('e0', calls=['e1'])])
SRV = task('t1', 'fcfs', [entry('e1')])

OR_BAD = '''      <task-activities>
        <activity name="a1" bound-to-entry="e1" host-demand-mean="1"/>
        <activity name="a2" host-demand-mean="1"/>
        <activity name="a3" host-demand-mean="1"/>
        <precedence>
          <pre><activity name="a1"/></pre>
          <post-OR><activity name="a2" prob="%s"/><activity name="a3" prob="%s"/></post-OR>
        </precedence>
        <reply-entry name="e1"><reply-activity name="a2"/></reply-entry>
      </task-activities>
'''

CASES = [
    ('dup_processor', proc('p0', [REF], 'inf') + proc('p0', [SRV]),
     'Duplicate processor name "p0".'),
    ('dup_task', proc('p0', [REF], 'inf') + proc('p1', [task('t0', 'fcfs', [entry('e1')])]),
     'Duplicate task name "t0".'),
    ('dup_entry', proc('p0', [REF], 'inf') + proc('p1', [task('t1', 'fcfs', [entry('e0')])]),
     'Duplicate entry name "e0".'),
    ('dup_activity',
     proc('p0', [REF], 'inf') +
     proc('p1', [task('t1', 'fcfs', [entry('e1')],
                      '      <task-activities>\n'
                      '        <activity name="e1_ph1" host-demand-mean="1"/>\n'
                      '      </task-activities>\n')]),
     'Duplicate activity name "e1_ph1" in task "t1".'),
    ('no_entries', proc('p0', [REF], 'inf') + proc('p1', [SRV]) +
     proc('p2', [task('t2', 'fcfs', [])]),
     'Task "t2" has no entries.'),
    ('no_reference_task',
     proc('p0', [task('t0', 'fcfs', [entry('e0', calls=['e1'])])]) + proc('p1', [SRV]),
     'The model has no reference task and no open arrivals.'),
    ('ref_receiver',
     proc('p0', [task('t0', 'ref', [entry('e0', calls=['e1'])])], 'inf') +
     proc('p1', [task('t1', 'fcfs', [entry('e1', calls=['e0'])])]),
     'Entry "e0" belongs to reference task "t0" and cannot receive requests.'),
    ('ref_replies',
     proc('p0', [task('t0', 'ref', [entry('e0', calls=['e1'])],
                      '      <task-activities>\n'
                      '        <reply-entry name="e0">\n'
                      '          <reply-activity name="e0_ph1"/>\n'
                      '        </reply-entry>\n'
                      '      </task-activities>\n')], 'inf') + proc('p1', [SRV]),
     'Entry "e0" belongs to reference task "t0" and cannot be replied to.'),
    ('ref_forwarding',
     proc('p0', [task('t0', 'ref', [entry('e0', calls=['e1'], fwd=[('e1', '0.5')])])], 'inf') +
     proc('p1', [SRV]),
     'Entry "e0" belongs to reference task "t0" and cannot forward requests.'),
    ('ref_open_arrivals',
     proc('p0', [task('t0', 'ref', [entry('e0', calls=['e1'], arrival='0.5')])], 'inf') +
     proc('p1', [SRV]),
     'Entry "e0" belongs to reference task "t0" and cannot have open arrivals.'),
    ('forwarding_probability_negative',
     proc('p0', [REF], 'inf') +
     proc('p1', [task('t1', 'fcfs', [entry('e1', fwd=[('e2', '-0.5')])])]) +
     proc('p2', [task('t2', 'fcfs', [entry('e2')])]),
     'Forwarding from entry "e1" to entry "e2" has an invalid probability of -0.5.'),
    ('forwarding_probability_total',
     proc('p0', [REF], 'inf') +
     proc('p1', [task('t1', 'fcfs', [entry('e1', fwd=[('e2', '0.7'), ('e3', '0.7')])])]) +
     proc('p2', [task('t2', 'fcfs', [entry('e2'), entry('e3')])]),
     'Entry "e1" has a total forwarding probability of 1.4.'),
    ('or_branch_probability_invalid',
     proc('p0', [REF], 'inf') +
     proc('p1', [task('t1', 'fcfs', ['      <entry name="e1" type="NONE"/>\n'],
                      OR_BAD % ('1.4', '-0.4'))]),
     'Activity "a2" in task "t1" has an invalid branch probability of 1.4.'),
    ('or_branch_probabilities_sum',
     proc('p0', [REF], 'inf') +
     proc('p1', [task('t1', 'fcfs', ['      <entry name="e1" type="NONE"/>\n'],
                      OR_BAD % ('0.4', '0.4'))]),
     'Branch probabilities of an OR-fork in task "t1" sum to 0.8 instead of 1.'),
]


def write_model(tmp_path, name, body):
    path = tmp_path / ('%s.lqnx' % name)
    path.write_text(HEAD + body + TAIL)
    return str(path)


@pytest.mark.parametrize('name,body,message', CASES, ids=[c[0] for c in CASES])
def test_defective_document_is_refused(tmp_path, name, body, message):
    path = write_model(tmp_path, name, body)
    with pytest.raises(Exception) as excinfo:
        LayeredNetwork.parse_xml(path)
    assert message in str(excinfo.value)


def test_consistent_document_is_accepted(tmp_path):
    path = write_model(tmp_path, 'good', proc('p0', [REF], 'inf') + proc('p1', [SRV]))
    model = LayeredNetwork.parse_xml(path)
    assert len(model.tasks) == 2
