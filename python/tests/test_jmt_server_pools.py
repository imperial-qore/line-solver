"""
Regression test: the JMT export must carry the server pools and job parallelism.

The native Python JSIM writer emitted no serverNames/serversPerServerType/
serverCompatibilities/schedulingPolicy at all, so a heterogeneous model silently
simulated as one homogeneous pool. Neither did any codebase emit classParallelism
(JMT's Server.serverNumRequired), the servers a job seizes for the whole of its
service.

JMT's SimLoader picks the Server constructor by the POSITIONAL types of the
parameters, so the five must follow ServiceStrategy and be emitted as one block:
written before it the engine dies with "Constructor of Section not found", and
there is no constructor taking classParallelism without the four pool parameters.
A station declaring parallelism alone therefore gets one synthetic pool carrying
all of its servers, since the pools, not maxJobs, size the pool once any exists.
"""
import xml.etree.ElementTree as ET

import pytest

from line_solver import (Exp, HeteroSchedPolicy, Network, OpenClass, Queue,
                         SchedStrategy, ServerType, Sink, SolverJMTOptions, Source)
from line_solver.api.solvers.jmt import handler


def _hetero():
    model = Network('hetero')
    source = Source(model, 'Source')
    queue = Queue(model, 'HeteroQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    jobs = OpenClass(model, 'Jobs')
    source.setArrival(jobs, Exp(3.0))
    queue.setService(jobs, Exp(1.5))
    queue.addServerType(ServerType('Fast', 2, [jobs]))
    queue.addServerType(ServerType('Slow', 3, [jobs]))
    queue.setHeteroSchedPolicy(HeteroSchedPolicy.FSF)
    model.link(Network.serialRouting(source, queue, sink))
    return model, queue, jobs


def _parallel(n, servers=4):
    model = Network('parallelism')
    source = Source(model, 'Source')
    queue = Queue(model, 'ParQueue', SchedStrategy.FCFS)
    sink = Sink(model, 'Sink')
    jobs = OpenClass(model, 'Jobs')
    source.setArrival(jobs, Exp(1.0))
    queue.setNumberOfServers(servers)
    queue.setService(jobs, Exp(10.0))
    if n > 1:
        queue.setServerParallelism(jobs, n)
    model.link(Network.serialRouting(source, queue, sink))
    return model, queue, jobs


def _server_section(model, tmp_path, name):
    path = str(tmp_path / (name + '.jsimg'))
    handler._write_jsim_file(model.getStruct(), path, SolverJMTOptions(), model)
    root = ET.parse(path).getroot()
    for section in root.iter('section'):
        if section.get('className') == 'Server':
            return section, open(path).read()
    raise AssertionError('no Server section was written')


def _names(section):
    return [p.get('name') for p in section.findall('parameter')]


def test_hetero_pools_are_exported_after_the_service_strategies(tmp_path):
    model, _, _ = _hetero()
    section, _ = _server_section(model, tmp_path, 'hetero')
    names = _names(section)
    for key in ('ServiceStrategy', 'classParallelism', 'serverNames',
                'serversPerServerType', 'serverCompatibilities', 'schedulingPolicy'):
        assert key in names, '%s must be exported' % key
    order = [names.index(k) for k in ('ServiceStrategy', 'classParallelism', 'serverNames',
                                      'serversPerServerType', 'serverCompatibilities',
                                      'schedulingPolicy')]
    assert order == sorted(order), 'the pool block must follow ServiceStrategy in this order'

    counts = section.find("parameter[@name='serversPerServerType']")
    assert [v.text for v in counts.iter('value')] == ['2', '3']
    pools = section.find("parameter[@name='serverNames']")
    assert [v.text for v in pools.iter('value')] == ['Fast', 'Slow']
    policy = section.find("parameter[@name='schedulingPolicy']")
    assert policy.find('value').text == 'FSF (Fastest Servers First)'
    # The pools size the station, so maxJobs is their sum
    assert section.find("parameter[@name='maxJobs']").find('value').text == '5'


def test_parallelism_alone_gets_a_synthetic_pool(tmp_path):
    model, _, _ = _parallel(2)
    section, _ = _server_section(model, tmp_path, 'par2')
    names = _names(section)
    assert 'classParallelism' in names
    par = section.find("parameter[@name='classParallelism']")
    assert [v.text for v in par.iter('value')] == ['2']
    # One synthetic pool, holding every server of the station
    counts = section.find("parameter[@name='serversPerServerType']")
    assert [v.text for v in counts.iter('value')] == ['4']
    pool_names = section.find("parameter[@name='serverNames']")
    assert [v.text for v in pool_names.iter('value')] == ['ParQueue - Server Type 1']


def test_a_plain_station_writes_no_pool_block(tmp_path):
    model, _, _ = _parallel(1)
    section, _ = _server_section(model, tmp_path, 'par1')
    names = _names(section)
    assert 'classParallelism' not in names
    assert 'serverNames' not in names
    assert 'schedulingPolicy' not in names


def test_parallelism_is_declared_and_bounded():
    model, queue, jobs = _parallel(1, servers=2)
    assert queue.getServerParallelism(jobs) == 1
    assert not queue.hasServerParallelism()
    assert not model.getUsedLangFeatures().list['ServerParallelism']

    queue.setServerParallelism(jobs, 2)
    assert queue.hasServerParallelism()
    assert model.getUsedLangFeatures().list['ServerParallelism']

    # A job seizing more servers than the station has could never enter service
    with pytest.raises(ValueError):
        queue.setServerParallelism(jobs, 3)
    with pytest.raises(ValueError):
        queue.setServerParallelism(jobs, 0)


def test_parallelism_survives_the_json_round_trip(tmp_path):
    from line_solver.io.linemodel_io import save_model, load_model

    model, queue, jobs = _parallel(3)
    path = str(tmp_path / 'par.json')
    save_model(model, path)
    back = load_model(path)
    q2 = back.getNodeByName('ParQueue')
    c2 = back.getClassByName('Jobs')
    assert q2.getServerParallelism(c2) == 3
