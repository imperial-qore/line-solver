"""
Workflow pattern detection and analysis functions.

This module provides functions for analyzing workflow patterns in
process models and queueing networks.
"""

import jpype
from .. import jlineMatrixFromArray, jlineMatrixToArray


def wf_analyzer(workflow_model):
    """
    Analyzes workflow structure and patterns.

    Args:
        workflow_model: Workflow model representation

    Returns:
        dict: Analysis results including patterns and metrics
    """
    result = jpype.JPackage('jline').api.wf.Wf_analyzerKt.wf_analyzer(
        jlineMatrixFromArray(workflow_model)
    )
    return {
        'patterns': jlineMatrixToArray(result.getPatterns()),
        'metrics': jlineMatrixToArray(result.getMetrics()),
        'structure': result.getStructure() if hasattr(result, 'getStructure') else None
    }


def wf_auto_integration(workflow_models):
    """
    Automatically integrates multiple workflow models.

    Args:
        workflow_models: List of workflow models to integrate

    Returns:
        dict: Integrated workflow model
    """
    result = jpype.JPackage('jline').api.wf.Wf_auto_integrationKt.wf_auto_integration(
        jlineMatrixFromArray(workflow_models)
    )
    return {
        'integrated_model': jlineMatrixToArray(result.getIntegratedModel()),
        'integration_points': jlineMatrixToArray(result.getIntegrationPoints()),
        'conflicts': result.getConflicts() if hasattr(result, 'getConflicts') else []
    }


def wf_branch_detector(workflow_model):
    """
    Detects branching patterns in workflow.

    Args:
        workflow_model: Workflow model representation

    Returns:
        dict: Branch detection results
    """
    result = jpype.JPackage('jline').api.wf.Wf_branch_detectorKt.wf_branch_detector(
        jlineMatrixFromArray(workflow_model)
    )
    return {
        'branches': jlineMatrixToArray(result.getBranches()),
        'branch_types': result.getBranchTypes() if hasattr(result, 'getBranchTypes') else [],
        'branch_probabilities': jlineMatrixToArray(result.getBranchProbabilities()) if hasattr(result, 'getBranchProbabilities') else None
    }


def wf_loop_detector(workflow_model):
    """
    Detects loop patterns in workflow.

    Args:
        workflow_model: Workflow model representation

    Returns:
        dict: Loop detection results
    """
    result = jpype.JPackage('jline').api.wf.Wf_loop_detectorKt.wf_loop_detector(
        jlineMatrixFromArray(workflow_model)
    )
    return {
        'loops': jlineMatrixToArray(result.getLoops()),
        'loop_bounds': jlineMatrixToArray(result.getLoopBounds()) if hasattr(result, 'getLoopBounds') else None,
        'nested_loops': result.getNestedLoops() if hasattr(result, 'getNestedLoops') else []
    }


def wf_parallel_detector(workflow_model):
    """
    Detects parallel execution patterns in workflow.

    Args:
        workflow_model: Workflow model representation

    Returns:
        dict: Parallel pattern detection results
    """
    result = jpype.JPackage('jline').api.wf.Wf_parallel_detectorKt.wf_parallel_detector(
        jlineMatrixFromArray(workflow_model)
    )
    return {
        'parallel_blocks': jlineMatrixToArray(result.getParallelBlocks()),
        'fork_points': jlineMatrixToArray(result.getForkPoints()),
        'join_points': jlineMatrixToArray(result.getJoinPoints()),
        'max_parallelism': int(result.getMaxParallelism()) if hasattr(result, 'getMaxParallelism') else None
    }


def wf_pattern_updater(workflow_model, patterns):
    """
    Updates workflow model with new patterns.

    Args:
        workflow_model: Original workflow model
        patterns: New patterns to apply

    Returns:
        dict: Updated workflow model
    """
    result = jpype.JPackage('jline').api.wf.Wf_pattern_updaterKt.wf_pattern_updater(
        jlineMatrixFromArray(workflow_model), jlineMatrixFromArray(patterns)
    )
    return {
        'updated_model': jlineMatrixToArray(result.getUpdatedModel()),
        'changes': result.getChanges() if hasattr(result, 'getChanges') else [],
        'validation': result.getValidation() if hasattr(result, 'getValidation') else None
    }


def wf_sequence_detector(workflow_model):
    """
    Detects sequential execution patterns in workflow.

    Args:
        workflow_model: Workflow model representation

    Returns:
        dict: Sequence detection results
    """
    result = jpype.JPackage('jline').api.wf.Wf_sequence_detectorKt.wf_sequence_detector(
        jlineMatrixFromArray(workflow_model)
    )
    return {
        'sequences': jlineMatrixToArray(result.getSequences()),
        'sequence_lengths': jlineMatrixToArray(result.getSequenceLengths()),
        'critical_path': jlineMatrixToArray(result.getCriticalPath()) if hasattr(result, 'getCriticalPath') else None
    }