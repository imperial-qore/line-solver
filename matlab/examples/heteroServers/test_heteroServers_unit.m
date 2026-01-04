%% Unit Tests for Heterogeneous Server Support
%
% This script contains unit tests for:
% - ServerType class
% - HeteroSchedPolicy enum
% - Queue heterogeneous server methods
% - NetworkStruct heterogeneous fields
%
% Run with: run('examples/heteroServers/test_heteroServers_unit.m')

clear;
fprintf('=== Unit Tests for Heterogeneous Servers ===\n\n');

testsPassed = 0;
testsFailed = 0;

%% Test 1: ServerType Creation
fprintf('Test 1: ServerType Creation... ');
try
    st = ServerType('Fast', 2);
    assert(strcmp(st.getName(), 'Fast'), 'Name should be Fast');
    assert(st.numOfServers == 2, 'Should have 2 servers');
    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 2: ServerType Set Number of Servers
fprintf('Test 2: ServerType Set Number of Servers... ');
try
    st = ServerType('Test', 1);
    assert(st.numOfServers == 1, 'Initial count should be 1');
    st.numOfServers = 5;
    assert(st.numOfServers == 5, 'Count should be 5');
    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 3: ServerType Compatibility
fprintf('Test 3: ServerType Compatibility... ');
try
    model = Network('TestModel');
    source = Source(model, 'Source');
    sink = Sink(model, 'Sink');

    class1 = OpenClass(model, 'Class1');
    class2 = OpenClass(model, 'Class2');

    st = ServerType('Fast', 2);

    % Initially not compatible
    assert(~st.isCompatible(class1), 'Should not be compatible with class1');
    assert(~st.isCompatible(class2), 'Should not be compatible with class2');

    % Add class1
    st.addCompatible(class1);
    assert(st.isCompatible(class1), 'Should be compatible with class1');
    assert(~st.isCompatible(class2), 'Should not be compatible with class2');

    % Add class2
    st.addCompatible(class2);
    assert(st.isCompatible(class1), 'Should be compatible with class1');
    assert(st.isCompatible(class2), 'Should be compatible with class2');

    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 4: ServerType setCompatible with Multiple Classes
fprintf('Test 4: ServerType setCompatible with Multiple Classes... ');
try
    model = Network('TestModel');
    source = Source(model, 'Source');
    sink = Sink(model, 'Sink');

    class1 = OpenClass(model, 'Class1');
    class2 = OpenClass(model, 'Class2');

    st = ServerType('Fast', 2);
    st.setCompatible({class1, class2});

    assert(st.isCompatible(class1), 'Should be compatible with class1');
    assert(st.isCompatible(class2), 'Should be compatible with class2');

    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 5: HeteroSchedPolicy Values
fprintf('Test 5: HeteroSchedPolicy Values... ');
try
    assert(HeteroSchedPolicy.ORDER == 0, 'ORDER should be 0');
    assert(HeteroSchedPolicy.ALIS == 1, 'ALIS should be 1');
    assert(HeteroSchedPolicy.ALFS == 2, 'ALFS should be 2');
    assert(HeteroSchedPolicy.FAIRNESS == 3, 'FAIRNESS should be 3');
    assert(HeteroSchedPolicy.FSF == 4, 'FSF should be 4');
    assert(HeteroSchedPolicy.RAIS == 5, 'RAIS should be 5');
    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 6: HeteroSchedPolicy toText
fprintf('Test 6: HeteroSchedPolicy toText... ');
try
    assert(strcmp(HeteroSchedPolicy.toText(HeteroSchedPolicy.FSF), 'FSF'), 'Should return FSF');
    assert(strcmp(HeteroSchedPolicy.toText(HeteroSchedPolicy.ALIS), 'ALIS'), 'Should return ALIS');
    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 7: Queue Initially Not Heterogeneous
fprintf('Test 7: Queue Initially Not Heterogeneous... ');
try
    model = Network('TestModel');
    queue = Queue(model, 'Queue', SchedStrategy.FCFS);

    assert(~queue.isHeterogeneous(), 'Queue should not be heterogeneous');
    assert(isempty(queue.serverTypes), 'Server types should be empty');

    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 8: Queue Add Server Type
fprintf('Test 8: Queue Add Server Type... ');
try
    model = Network('TestModel');
    queue = Queue(model, 'Queue', SchedStrategy.FCFS);

    fast = ServerType('Fast', 2);
    slow = ServerType('Slow', 3);

    queue.addServerType(fast);
    assert(queue.isHeterogeneous(), 'Queue should be heterogeneous');
    assert(length(queue.serverTypes) == 1, 'Should have 1 server type');

    queue.addServerType(slow);
    assert(length(queue.serverTypes) == 2, 'Should have 2 server types');

    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 9: Queue Set Hetero Sched Policy
fprintf('Test 9: Queue Set Hetero Sched Policy... ');
try
    model = Network('TestModel');
    queue = Queue(model, 'Queue', SchedStrategy.FCFS);

    queue.setHeteroSchedPolicy(HeteroSchedPolicy.FSF);
    assert(queue.heteroSchedPolicy == HeteroSchedPolicy.FSF, 'Policy should be FSF');

    queue.setHeteroSchedPolicy(HeteroSchedPolicy.ALIS);
    assert(queue.heteroSchedPolicy == HeteroSchedPolicy.ALIS, 'Policy should be ALIS');

    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 10: Queue Set Hetero Service
fprintf('Test 10: Queue Set Hetero Service... ');
try
    model = Network('TestModel');
    source = Source(model, 'Source');
    queue = Queue(model, 'Queue', SchedStrategy.FCFS);
    sink = Sink(model, 'Sink');

    jobClass = OpenClass(model, 'Jobs');
    source.setArrival(jobClass, Exp(1.0));

    fast = ServerType('Fast', 2);
    slow = ServerType('Slow', 3);
    fast.setCompatible(jobClass);
    slow.setCompatible(jobClass);

    queue.addServerType(fast);
    queue.addServerType(slow);

    queue.setHeteroService(jobClass, fast, Exp(2.0));
    queue.setHeteroService(jobClass, slow, Exp(1.0));

    assert(queue.isHeterogeneous(), 'Queue should be heterogeneous');
    assert(queue.heteroServiceDistributions.Count > 0, 'Should have service distributions');

    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 11: NetworkStruct Heterogeneous Fields
fprintf('Test 11: NetworkStruct Heterogeneous Fields... ');
try
    model = Network('TestModel');
    source = Source(model, 'Source');
    queue = Queue(model, 'Queue', SchedStrategy.FCFS);
    sink = Sink(model, 'Sink');

    jobClass = OpenClass(model, 'Jobs');
    source.setArrival(jobClass, Exp(1.0));

    fast = ServerType('Fast', 2);
    slow = ServerType('Slow', 3);
    fast.setCompatible(jobClass);
    slow.setCompatible(jobClass);

    queue.addServerType(fast);
    queue.addServerType(slow);
    queue.setHeteroSchedPolicy(HeteroSchedPolicy.FSF);
    queue.setHeteroService(jobClass, fast, Exp(2.0));
    queue.setHeteroService(jobClass, slow, Exp(1.0));

    model.link(Network.serialRouting({source, queue, sink}));

    sn = model.getStruct();
    ist = sn.nodeToStation(queue.index);

    assert(sn.nservertypes(ist) == 2, 'Should have 2 server types');
    assert(strcmp(sn.servertypenames{ist}{1}, 'Fast'), 'First type should be Fast');
    assert(strcmp(sn.servertypenames{ist}{2}, 'Slow'), 'Second type should be Slow');
    assert(sn.serverspertype{ist}(1) == 2, 'Fast should have 2 servers');
    assert(sn.serverspertype{ist}(2) == 3, 'Slow should have 3 servers');
    assert(sn.heteroschedpolicy(ist) == HeteroSchedPolicy.FSF, 'Policy should be FSF');

    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 12: Compatibility Matrix in NetworkStruct
fprintf('Test 12: Compatibility Matrix in NetworkStruct... ');
try
    model = Network('TestModel');
    source = Source(model, 'Source');
    queue = Queue(model, 'Queue', SchedStrategy.FCFS);
    sink = Sink(model, 'Sink');

    classA = OpenClass(model, 'ClassA');
    classB = OpenClass(model, 'ClassB');

    source.setArrival(classA, Exp(1.0));
    source.setArrival(classB, Exp(0.5));

    % Fast serves only ClassA
    fast = ServerType('Fast', 2);
    fast.setCompatible(classA);

    % Slow serves both
    slow = ServerType('Slow', 3);
    slow.setCompatible({classA, classB});

    queue.addServerType(fast);
    queue.addServerType(slow);
    queue.setHeteroService(classA, fast, Exp(2.0));
    queue.setHeteroService(classA, slow, Exp(1.0));
    queue.setHeteroService(classB, slow, Exp(1.5));

    P = model.initRoutingMatrix();
    P{classA} = Network.serialRouting({source, queue, sink});
    P{classB} = Network.serialRouting({source, queue, sink});
    model.link(P);

    sn = model.getStruct();
    ist = sn.nodeToStation(queue.index);

    compat = sn.servercompat{ist};
    % compat(type, class)
    % Fast (type 1): compatible with ClassA (1), not ClassB (0)
    % Slow (type 2): compatible with both (1, 1)
    assert(compat(1, 1) == 1, 'Fast should be compatible with ClassA');
    assert(compat(1, 2) == 0, 'Fast should not be compatible with ClassB');
    assert(compat(2, 1) == 1, 'Slow should be compatible with ClassA');
    assert(compat(2, 2) == 1, 'Slow should be compatible with ClassB');

    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Test 13: Total Server Count
fprintf('Test 13: Total Server Count... ');
try
    model = Network('TestModel');
    queue = Queue(model, 'Queue', SchedStrategy.FCFS);

    st1 = ServerType('Type1', 2);
    st2 = ServerType('Type2', 3);
    st3 = ServerType('Type3', 5);

    queue.addServerType(st1);
    queue.addServerType(st2);
    queue.addServerType(st3);

    totalServers = 0;
    for i = 1:length(queue.serverTypes)
        totalServers = totalServers + queue.serverTypes{i}.numOfServers;
    end

    assert(totalServers == 10, 'Total servers should be 10');

    fprintf('PASSED\n');
    testsPassed = testsPassed + 1;
catch ME
    fprintf('FAILED: %s\n', ME.message);
    testsFailed = testsFailed + 1;
end

%% Summary
fprintf('\n=== Test Summary ===\n');
fprintf('Passed: %d\n', testsPassed);
fprintf('Failed: %d\n', testsFailed);
fprintf('Total:  %d\n', testsPassed + testsFailed);

if testsFailed == 0
    fprintf('\nAll tests passed!\n');
else
    fprintf('\nSome tests failed. Please review.\n');
end
