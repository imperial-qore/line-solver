"""
Load WfCommons JSON workflows into LINE Workflow objects.

WfCommons (https://github.com/wfcommons/workflow-schema) is a standard format
for scientific workflow traces. ``WfCommonsLoader`` reads a trace into a
:class:`line_solver.lang.workflow.Workflow` for queueing analysis: each task
becomes an activity whose host demand is fitted to its recorded runtime, an
AND-fork/AND-join pair is detected wherever every child of a task feeds one
common task directly, and every remaining edge becomes a serial precedence.

Supported schema versions: 1.3, 1.4, 1.5. In schema 1.4 and earlier the tasks
sit under ``workflow.tasks`` with the runtime embedded in each task; from 1.5
they sit under ``workflow.specification.tasks`` with the runtimes in a separate
``workflow.execution.tasks`` section.

Example:
    >>> wf = WfCommonsLoader.load('workflow.json')
    >>> alpha, T = wf.toPH()

    >>> wf = WfCommonsLoader.load('workflow.json',
    ...                           {'distributionType': 'exp', 'defaultRuntime': 1.0})

Port of:
    matlab/src/io/WfCommonsLoader.m (reference)
    jar/src/main/java/jline/io/WfCommonsLoader.java, WfCommonsOptions.java

Copyright (c) 2012-2026, Imperial College London
All rights reserved.
"""

import json
import os
import re
import urllib.request
from collections import deque
from enum import Enum
from typing import Any, Dict, List, Mapping, Optional, Union

from ..constants import GlobalConstants
from ..api.io.logging import line_error, line_warning


class WfCommonsOptions:
    """
    Options of :class:`WfCommonsLoader`, the twin of Java ``WfCommonsOptions``
    and of the MATLAB options struct.

    Attributes (MATLAB field names):
        distributionType: 'exp' (default), 'det', 'aph' or 'hyperexp'
        defaultSCV: SCV used by 'aph' and 'hyperexp' (default 1.0)
        defaultRuntime: runtime of a task with no recorded one (default 1.0)
        useExecutionData: read the recorded runtimes when present (default True)
        storeMetadata: attach the WfCommons task metadata to each activity (default True)
        resolveChildrenByName: look a child reference up among the task ids first and,
            failing that, among the task names (default True). Schema 1.4 Pegasus traces
            (e.g. Montage) give each task an id such as 'ID0000001' but list its children
            by name, so an id-only lookup drops every edge. A name shared by several tasks
            is ambiguous and does not resolve; a reference resolving to no task is reported
            by line_warning. False gives the id-only lookup.

    Edges are the union of the ``children`` and ``parents`` lists, each (from, to) pair once and
    the children-derived edges first, so a trace that gives either form (the Montage dss trace
    lists parents only) yields the same graph. Parent references resolve by the same rule.

    The setters return the object, so options chain as in Java.
    """

    class DistributionType(Enum):
        EXP = 'exp'
        DET = 'det'
        APH = 'aph'
        HYPEREXP = 'hyperexp'

    FIELDS = ('distributionType', 'defaultSCV', 'defaultRuntime',
              'useExecutionData', 'storeMetadata', 'resolveChildrenByName')

    def __init__(self, distributionType: Union[str, 'WfCommonsOptions.DistributionType'] = 'exp',
                 defaultSCV: float = 1.0, defaultRuntime: float = 1.0,
                 useExecutionData: bool = True, storeMetadata: bool = True,
                 resolveChildrenByName: bool = True):
        self.distributionType = distributionType
        self.defaultSCV = defaultSCV
        self.defaultRuntime = defaultRuntime
        self.useExecutionData = useExecutionData
        self.storeMetadata = storeMetadata
        self.resolveChildrenByName = resolveChildrenByName

    @property
    def distributionType(self) -> str:
        return self._distributionType

    @distributionType.setter
    def distributionType(self, value):
        if isinstance(value, WfCommonsOptions.DistributionType):
            value = value.value
        self._distributionType = str(value)

    # Java-style accessors
    def getDistributionType(self) -> str:
        return self.distributionType

    def setDistributionType(self, value) -> 'WfCommonsOptions':
        self.distributionType = value
        return self

    def getDefaultSCV(self) -> float:
        return self.defaultSCV

    def setDefaultSCV(self, value: float) -> 'WfCommonsOptions':
        self.defaultSCV = value
        return self

    def getDefaultRuntime(self) -> float:
        return self.defaultRuntime

    def setDefaultRuntime(self, value: float) -> 'WfCommonsOptions':
        self.defaultRuntime = value
        return self

    def isUseExecutionData(self) -> bool:
        return self.useExecutionData

    def setUseExecutionData(self, value: bool) -> 'WfCommonsOptions':
        self.useExecutionData = value
        return self

    def isStoreMetadata(self) -> bool:
        return self.storeMetadata

    def setStoreMetadata(self, value: bool) -> 'WfCommonsOptions':
        self.storeMetadata = value
        return self

    def isResolveChildrenByName(self) -> bool:
        return self.resolveChildrenByName

    def setResolveChildrenByName(self, value: bool) -> 'WfCommonsOptions':
        self.resolveChildrenByName = value
        return self

    @staticmethod
    def exponential() -> 'WfCommonsOptions':
        """Options fitting every runtime with an exponential law."""
        return WfCommonsOptions(distributionType='exp')

    @staticmethod
    def deterministic() -> 'WfCommonsOptions':
        """Options fitting every runtime with a deterministic law."""
        return WfCommonsOptions(distributionType='det')

    def to_dict(self) -> Dict[str, Any]:
        """The options as a MATLAB-style struct (dict keyed by MATLAB field name)."""
        return {f: getattr(self, f) for f in WfCommonsOptions.FIELDS}

    def __repr__(self) -> str:
        return 'WfCommonsOptions(%s)' % ', '.join('%s=%r' % kv for kv in self.to_dict().items())


class WfCommonsLoader:
    """
    Static loader of WfCommons JSON workflows, the port of MATLAB
    ``WfCommonsLoader``. Every method takes an optional ``options`` argument,
    given as a :class:`WfCommonsOptions`, as a dict with the MATLAB field names
    (a partial dict takes the defaults for the missing fields), or None.
    """

    SUPPORTED_SCHEMA_VERSIONS = ('1.3', '1.4', '1.5')

    # ------------------------------------------------------------------ public

    @staticmethod
    def load(jsonFile: str, options: Optional[Union[WfCommonsOptions, Mapping[str, Any]]] = None):
        """
        Load a WfCommons JSON file into a Workflow.

        Args:
            jsonFile: path to the WfCommons JSON file
            options: loader options (see :class:`WfCommonsOptions`)

        Returns:
            Workflow
        """
        options = WfCommonsLoader.parseOptions(options)
        with open(jsonFile, 'r') as f:
            data = json.load(f)
        WfCommonsLoader._validate_schema(data)
        name = WfCommonsLoader._extract_name(data, jsonFile)
        return WfCommonsLoader._build_workflow(data, name, options)

    @staticmethod
    def loadFromStruct(data: Mapping[str, Any],
                       options: Optional[Union[WfCommonsOptions, Mapping[str, Any]]] = None):
        """
        Load a Workflow from an already parsed WfCommons document (the dict that
        ``json.load`` returns, the counterpart of MATLAB's jsondecode struct).
        """
        options = WfCommonsLoader.parseOptions(options)
        WfCommonsLoader._validate_schema(data)
        name = WfCommonsLoader._extract_name(data, 'struct_input')
        return WfCommonsLoader._build_workflow(data, name, options)

    @staticmethod
    def loadFromString(json_str: str,
                       options: Optional[Union[WfCommonsOptions, Mapping[str, Any]]] = None):
        """Load a Workflow from a WfCommons JSON string (Java ``loadFromString``)."""
        options = WfCommonsLoader.parseOptions(options)
        data = json.loads(json_str)
        WfCommonsLoader._validate_schema(data)
        name = WfCommonsLoader._extract_name(data, 'string_input')
        return WfCommonsLoader._build_workflow(data, name, options)

    @staticmethod
    def loadFromUrl(urlString: str,
                    options: Optional[Union[WfCommonsOptions, Mapping[str, Any]]] = None):
        """
        Load a WfCommons JSON document from a URL, e.g. a file of the
        wfcommons/pegasus-instances repository. The workflow is named after the
        document, or after the last path component of the URL when it has none.
        """
        options = WfCommonsLoader.parseOptions(options)
        with urllib.request.urlopen(urlString) as resp:
            data = json.loads(resp.read().decode('utf-8'))
        WfCommonsLoader._validate_schema(data)
        default_name = os.path.splitext(urlString.rstrip('/').rsplit('/', 1)[-1])[0]
        name = WfCommonsLoader._extract_name(data, default_name)
        return WfCommonsLoader._build_workflow(data, name, options)

    @staticmethod
    def validateFile(jsonFile: str) -> bool:
        """True when the file parses as JSON and passes the WfCommons schema checks."""
        try:
            with open(jsonFile, 'r') as f:
                data = json.load(f)
            WfCommonsLoader._validate_schema(data)
            return True
        except Exception:
            return False

    @staticmethod
    def defaultOptions() -> WfCommonsOptions:
        """The default options, as MATLAB parseOptions fills them."""
        return WfCommonsOptions()

    @staticmethod
    def parseOptions(options: Optional[Union[WfCommonsOptions, Mapping[str, Any]]] = None) -> WfCommonsOptions:
        """
        Normalise ``options`` to a :class:`WfCommonsOptions`, filling the missing
        fields with their defaults as MATLAB ``parseOptions`` does.
        """
        if options is None:
            return WfCommonsOptions()
        if isinstance(options, WfCommonsOptions):
            return options
        if isinstance(options, Mapping):
            unknown = [k for k in options if k not in WfCommonsOptions.FIELDS]
            if unknown:
                line_warning('WfCommonsLoader', 'Ignoring unknown option(s): %s.', ', '.join(map(str, unknown)))
            return WfCommonsOptions(**{k: v for k, v in options.items() if k in WfCommonsOptions.FIELDS})
        raise TypeError('WfCommonsLoader options must be a WfCommonsOptions, a dict or None, not %s.'
                        % type(options).__name__)

    # snake_case aliases
    load_from_struct = loadFromStruct
    load_from_string = loadFromString
    load_from_url = loadFromUrl
    validate_file = validateFile
    default_options = defaultOptions
    parse_options = parseOptions

    # ----------------------------------------------------------------- private

    @staticmethod
    def _validate_schema(data: Any) -> None:
        if not isinstance(data, Mapping) or 'schemaVersion' not in data:
            line_error('WfCommonsLoader', 'Missing schemaVersion field in WfCommons JSON.')
        version = str(data['schemaVersion'])
        if version not in WfCommonsLoader.SUPPORTED_SCHEMA_VERSIONS:
            line_warning('WfCommonsLoader', 'Schema version %s may not be fully supported.', version)
        if 'workflow' not in data:
            line_error('WfCommonsLoader', 'Missing workflow field in WfCommons JSON.')
        # Schema 1.4 and earlier use workflow.tasks, schema 1.5+ workflow.specification.tasks
        workflow = data['workflow'] or {}
        if 'specification' in workflow:
            tasks = (workflow['specification'] or {}).get('tasks')
        else:
            tasks = workflow.get('tasks')
        if not tasks:
            line_error('WfCommonsLoader', 'Workflow must have at least one task.')

    @staticmethod
    def _extract_name(data: Mapping[str, Any], default_name: str) -> str:
        name = data.get('name')
        if not name:
            name = os.path.splitext(os.path.basename(default_name.replace('\\', '/')))[0]
        name = re.sub(r'[^a-zA-Z0-9_]', '_', str(name))
        return name if name else 'Workflow'

    @staticmethod
    def _task_id(task: Mapping[str, Any], idx: int) -> str:
        # Schema 1.4 may identify a task by 'name' only; idx is one-based as in MATLAB
        if 'id' in task:
            return task['id']
        if 'name' in task:
            return task['name']
        return 'task_%d' % idx

    @staticmethod
    def _as_list(value: Any) -> List[Any]:
        if value is None:
            return []
        if isinstance(value, str):
            return [value]
        return list(value)

    @staticmethod
    def _build_workflow(data: Mapping[str, Any], workflow_name: str, options: WfCommonsOptions):
        from ..lang.workflow import Workflow

        wf = Workflow(workflow_name)
        workflow = data['workflow']
        is_legacy = 'specification' not in workflow
        tasks = list(workflow['tasks'] if is_legacy else workflow['specification']['tasks'])
        task_id = WfCommonsLoader._task_id

        # Execution data, keyed by task id
        exec_map: Dict[str, Mapping[str, Any]] = {}
        if options.useExecutionData:
            if is_legacy:
                # Schema 1.4: the runtime is embedded in the task objects
                for i, task in enumerate(tasks):
                    exec_map[task_id(task, i + 1)] = task
            elif 'execution' in workflow and 'tasks' in (workflow['execution'] or {}):
                for et in WfCommonsLoader._as_list(workflow['execution']['tasks']):
                    exec_map[et['id']] = et

        # Phase 1: activities
        task_map = {}
        task_idx = {}
        for i, task in enumerate(tasks):
            tid = task_id(task, i + 1)
            runtime = options.defaultRuntime
            ed = exec_map.get(tid)
            if ed is not None and 'runtimeInSeconds' in ed:
                runtime = ed['runtimeInSeconds']
            act = wf.addActivity(tid, WfCommonsLoader._fit_distribution(float(runtime), options))
            task_map[tid] = act
            task_idx[tid] = i
            if options.storeMetadata:
                act.metadata = WfCommonsLoader._extract_metadata(task, tid, exec_map)

        # Phase 2: adjacency (zero-based task indices); a child or parent resolves by task id, then by a unique task name.
        # Edges are the union of `children` (task -> child) and `parents` (parent -> task), each pair once, so a trace
        # listing either form (the Montage dss trace lists parents only) gives the same graph.
        n = len(tasks)
        by_name = bool(options.resolveChildrenByName)
        name_idx: Dict[str, int] = {}
        name_count: Dict[str, int] = {}
        if by_name:
            for i, task in enumerate(tasks):
                nm = task.get('name')
                if isinstance(nm, str):
                    name_idx[nm] = i
                    name_count[nm] = name_count.get(nm, 0) + 1
        adj: List[List[int]] = [[] for _ in range(n)]
        in_deg = [0] * n
        out_deg = [0] * n
        unresolved: List[str] = []
        seen = set()
        for field in ('children', 'parents'):
            for i, task in enumerate(tasks):
                for ref in WfCommonsLoader._as_list(task.get(field)):
                    ref = str(ref)
                    if ref in task_idx:
                        r = task_idx[ref]
                    elif by_name and name_count.get(ref) == 1:
                        r = name_idx[ref]
                    else:
                        unresolved.append(ref)
                        continue
                    edge = (i, r) if field == 'children' else (r, i)
                    if edge in seen:
                        continue
                    seen.add(edge)
                    adj[edge[0]].append(edge[1])
                    out_deg[edge[0]] += 1
                    in_deg[edge[1]] += 1
        if unresolved:
            shown = list(dict.fromkeys(unresolved))
            line_warning('WfCommonsLoader', '%d child/parent reference(s) match no task id%s and were dropped: %s%s',
                         len(unresolved), ' or unique task name' if by_name else '',
                         ', '.join(shown[:5]), ', ...' if len(shown) > 5 else '')

        # Phase 3: precedences
        acts = [task_map[task_id(t, i + 1)] for i, t in enumerate(tasks)]
        WfCommonsLoader._add_precedences(wf, acts, adj, in_deg, out_deg)
        return wf

    @staticmethod
    def _fit_distribution(runtime: float, options: WfCommonsOptions):
        from ..distributions import Exp, Det, APH, HyperExp, Immediate

        if runtime <= GlobalConstants.FineTol:
            return Immediate()
        kind = options.distributionType.lower()
        if kind == 'det':
            return Det(runtime)
        if kind == 'aph':
            return APH.fitMeanAndSCV(runtime, options.defaultSCV)
        if kind == 'hyperexp' and options.defaultSCV > 1.0:
            return HyperExp.fitMeanAndSCV(runtime, options.defaultSCV)
        # 'exp', 'hyperexp' with SCV <= 1, and any other name fall back to Exp, as in MATLAB
        return Exp.fitMean(runtime)

    @staticmethod
    def _extract_metadata(task: Mapping[str, Any], tid: str,
                          exec_map: Mapping[str, Mapping[str, Any]]) -> Dict[str, Any]:
        metadata: Dict[str, Any] = {'taskId': tid}
        for f in ('name', 'inputFiles', 'outputFiles'):
            if f in task:
                metadata[f] = task[f]
        ed = exec_map.get(tid)
        if ed is not None:
            for f in ('executedAt', 'command', 'coreCount', 'avgCPU', 'readBytes',
                      'writtenBytes', 'memoryInBytes', 'energyInKWh', 'avgPowerInW',
                      'priority', 'machines'):
                if f in ed:
                    metadata[f] = ed[f]
        return metadata

    @staticmethod
    def _add_precedences(wf, acts, adj, in_deg, out_deg) -> None:
        from ..lang.workflow import Workflow

        n = len(acts)
        processed = set()
        # Step 1: AND-fork/AND-join pairs
        for i in range(n):
            if out_deg[i] > 1:
                children = adj[i]
                join = WfCommonsLoader._find_common_join(children, adj, in_deg, n)
                if join is not None:
                    post = [acts[c] for c in children]
                    wf.addPrecedence(Workflow.AndFork(acts[i], post))
                    wf.addPrecedence(Workflow.AndJoin(post, acts[join]))
                    for c in children:
                        processed.add((i, c))
                        processed.add((c, join))
        # Step 2: every remaining edge is serial
        for i in range(n):
            for j in adj[i]:
                if (i, j) not in processed:
                    wf.addPrecedence(Workflow.Serial(acts[i], acts[j]))

    @staticmethod
    def _find_common_join(children, adj, in_deg, n) -> Optional[int]:
        if not children:
            return None
        common = None
        for c in children:
            r = WfCommonsLoader._reachable(c, adj, n)
            common = r if common is None else common & r
        # MATLAB intersect returns the common nodes in ascending order
        for node in sorted(common):
            if in_deg[node] >= len(children) and all(node in adj[c] for c in children):
                return node
        return None

    @staticmethod
    def _reachable(start: int, adj, n) -> set:
        visited = {start}
        reach = set()
        queue = deque([start])
        while queue:
            cur = queue.popleft()
            for nxt in adj[cur]:
                if nxt not in visited:
                    visited.add(nxt)
                    queue.append(nxt)
                    reach.add(nxt)
        return reach


__all__ = ['WfCommonsLoader', 'WfCommonsOptions']
