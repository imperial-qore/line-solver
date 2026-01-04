"""
Workflow (WF) Analysis.

Native Python implementations for workflow pattern detection and analysis
in queueing networks. Provides comprehensive tools for identifying workflow
structures including sequences, parallel patterns, loops, and branches.

This module implements the AUTO workflow analysis algorithms from the MDN
toolbox, enabling intelligent solver selection based on detected patterns.

Key Classes:
    WorkflowManager: Main facade for workflow management and optimization
    WfAnalyzer: Workflow analyzer for pattern detection and optimization
    WfAutoIntegration: AUTO solver integration with workflow analysis

Pattern Detectors:
    detect_sequences: Detect sequential patterns in workflow
    detect_parallel: Detect parallel execution patterns
    detect_loops: Detect loop/iteration patterns
    detect_branches: Detect branching decision patterns

Data Classes:
    WorkflowAnalysis: Complete workflow analysis results
    DetectedPatterns: All detected patterns
    BranchPattern: Branch pattern with probabilities
    ServiceParameters: Phase-type service parameters
    UpdatedWorkflow: Optimized workflow structure

References:
    Original Kotlin: jar/src/main/kotlin/jline/api/wf/
    MDN toolbox AUTO solver algorithms
"""

# Sequence detector
from .sequence_detector import (
    detect_sequences,
    validate_sequence,
    get_sequence_stats,
)

# Branch detector
from .branch_detector import (
    BranchPattern,
    detect_branches,
    validate_branch_pattern,
    calculate_branch_diversity,
    get_branch_stats,
    find_most_probable_branch,
    find_least_probable_branch,
)

# Loop detector
from .loop_detector import (
    detect_loops,
    get_loop_probability,
    validate_loop_pattern,
    get_expected_loop_iterations,
    get_loop_stats,
)

# Parallel detector
from .parallel_detector import (
    detect_parallel,
    validate_parallel_pattern,
    get_parallel_stats,
)

# Pattern updater
from .pattern_updater import (
    ServiceParameters,
    UpdatedWorkflow,
    update_patterns,
    validate_updated_workflow,
    get_update_stats,
)

# Analyzer
from .analyzer import (
    WorkflowRepresentation,
    DetectedPatterns,
    WorkflowAnalysis,
    WfAnalyzer,
)

# Auto integration
from .auto_integration import (
    ExtendedSolverRecommendation,
    WfAutoIntegration,
)

# Workflow manager
from .workflow_manager import (
    WorkflowAnalysisResult,
    WorkflowManager,
    quick_analysis,
    get_optimal_solver,
)

__all__ = [
    # Sequence detector
    'detect_sequences',
    'validate_sequence',
    'get_sequence_stats',
    # Branch detector
    'BranchPattern',
    'detect_branches',
    'validate_branch_pattern',
    'calculate_branch_diversity',
    'get_branch_stats',
    'find_most_probable_branch',
    'find_least_probable_branch',
    # Loop detector
    'detect_loops',
    'get_loop_probability',
    'validate_loop_pattern',
    'get_expected_loop_iterations',
    'get_loop_stats',
    # Parallel detector
    'detect_parallel',
    'validate_parallel_pattern',
    'get_parallel_stats',
    # Pattern updater
    'ServiceParameters',
    'UpdatedWorkflow',
    'update_patterns',
    'validate_updated_workflow',
    'get_update_stats',
    # Analyzer
    'WorkflowRepresentation',
    'DetectedPatterns',
    'WorkflowAnalysis',
    'WfAnalyzer',
    # Auto integration
    'ExtendedSolverRecommendation',
    'WfAutoIntegration',
    # Workflow manager
    'WorkflowAnalysisResult',
    'WorkflowManager',
    'quick_analysis',
    'get_optimal_solver',
]
