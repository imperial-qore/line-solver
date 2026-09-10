"""The .lqnx "pri" discipline is preemptive priority resume, not FCFS.

lqns writes SCHEDULE_PPR as "pri" (LQIO::SCHEDULE::PPR in lqiolib labels.cpp),
and PPR preempts: a more urgent arrival takes the server and the displaced job
resumes its residual work (lqns/processor.cc builds a PR_FCFS_Server). LINE
holds that as FCFSPRPRIO. This reader used to fall through to FCFS, which is
the silent case: an unrecognized discipline became an ordinary queue and the
answer looked plausible. MATLAB and the JAR refused the token outright and C++
read it as the NON-preemptive HOL, so one file was three different models.

"pp" is the spelling of lqn-core.xsd, which is stale: it appears nowhere in the
lqns 6.2.31 sources, so no file lqns produces carries it. It is accepted on the
read side and never written. The twins are cpp/tests/test_lqnx_pri_sched.cpp,
jar LqnxPriSchedTest.java and MATLAB test_lqnx_pri_sched.m.
"""
import os
import tempfile

from line_solver import LayeredNetwork, SchedStrategy


def document(proc_sched, task_sched):
    """A reference task on an inf host calling one server task on host PS."""
    return (
        '<?xml version="1.0"?>\n'
        '<lqn-model name="pri" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">\n'
        '  <processor name="PC" scheduling="inf">\n'
        '    <task name="TC" scheduling="ref" multiplicity="2">\n'
        '      <entry name="EC" type="PH1PH2">\n'
        '        <entry-phase-activities>\n'
        '          <activity name="EC_ph1" phase="1" host-demand-mean="1">\n'
        '            <synch-call dest="ES" calls-mean="1"/>\n'
        '          </activity>\n'
        '        </entry-phase-activities>\n'
        '      </entry>\n'
        '    </task>\n'
        '  </processor>\n'
        '  <processor name="PS" scheduling="%s">\n'
        '    <task name="TS" scheduling="%s" priority="3">\n'
        '      <entry name="ES" type="PH1PH2">\n'
        '        <entry-phase-activities>\n'
        '          <activity name="ES_ph1" phase="1" host-demand-mean="2"/>\n'
        '        </entry-phase-activities>\n'
        '      </entry>\n'
        '    </task>\n'
        '  </processor>\n'
        '</lqn-model>\n' % (proc_sched, task_sched))


def write(dirname, name, body):
    path = os.path.join(dirname, name + '.lqnx')
    with open(path, 'w') as fd:
        fd.write(body)
    return path


def host_sched(model, name):
    for h in model.processors:
        if h.name == name:
            return h.getScheduling()
    raise AssertionError('no host ' + name)


def task_sched(model, name):
    for h in model.processors:
        for t in h.tasks:
            if t.name == name:
                return t.getScheduling()
    raise AssertionError('no task ' + name)


def test_pri_is_preemptive_priority_resume():
    d = tempfile.mkdtemp(prefix='lqnx_pri_')
    model = LayeredNetwork.parseXML(write(d, 'pri', document('pri', 'pri')))
    assert host_sched(model, 'PS') == SchedStrategy.FCFSPRPRIO
    assert task_sched(model, 'TS') == SchedStrategy.FCFSPRPRIO
    assert host_sched(model, 'PS') != SchedStrategy.FCFS

    # the stale schema spelling reads the same way
    stale = LayeredNetwork.parseXML(write(d, 'pp', document('pp', 'pp')))
    assert host_sched(stale, 'PS') == SchedStrategy.FCFSPRPRIO
    assert task_sched(stale, 'TS') == SchedStrategy.FCFSPRPRIO


def test_writer_emits_pri_the_spelling_lqns_reads_back():
    d = tempfile.mkdtemp(prefix='lqnx_pri_rt_')
    model = LayeredNetwork.parseXML(write(d, 'in', document('pri', 'pri')))
    out = os.path.join(d, 'out.lqnx')
    model.writeXML(out)

    txt = open(out).read()
    assert 'scheduling="pri"' in txt
    # fcfsprprio is LINE's own name for the discipline and is not valid LQN
    assert 'fcfsprprio' not in txt

    back = LayeredNetwork.parseXML(out)
    assert host_sched(back, 'PS') == SchedStrategy.FCFSPRPRIO
    assert task_sched(back, 'TS') == SchedStrategy.FCFSPRPRIO
