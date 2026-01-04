"""
LINE I/O and Model Transformation Module.

This module provides input/output operations and model transformation utilities
for the LINE solver, including:
- Logging and output utilities (line_warning, line_error, line_printf, line_debug)
- XML import/export (NetworkXMLIO)
- JMT import/export (jmt2line, jsim2line, jmva2line, qn2jsimg)
- JMT viewers (jsimg_view, jsimw_view, jmva_view)
- Code generation (qn2python, line2python, qn2java, lqn2java, line2java)
- Model format converters (qn2line, qn2lqn, lqn2qn, mapqn2renv)
- Model adaptation (ModelAdapter, chain aggregation, tagging, fork-join transforms)

Port from:
    - matlab/src/io/

Example usage:
    from line_solver.api.io import (
        line_warning, line_printf, line_verbosity, VerboseLevel,
        NetworkXMLIO, qn2python, qn2java,
        ModelAdapter, aggregate_chains, mmt, ht,
        jmt2line, qn2jsimg, jsimg_view
    )

    # Set verbosity
    line_verbosity(VerboseLevel.DEBUG)

    # Print messages
    line_printf("Processing model...\\n")
    line_warning("my_function", "This is a warning: %s", "details")

    # Export model to XML
    NetworkXMLIO.export_to_xml(model, 'model.xml')

    # Generate Python code from model
    code = qn2python(model, 'my_model')

    # Generate Java code from model
    java_code = qn2java(model, 'my_model')

    # Export to JMT format and view
    jsimg_file = qn2jsimg(model)
    jsimg_view(jsimg_file)

    # Aggregate chains
    chain_model, alpha, deagg_info = aggregate_chains(model)

    # Transform fork-join networks
    mmt_result = mmt(model)
"""

# Logging utilities
from .logging import (
    VerboseLevel,
    LineLogger,
    line_printf,
    line_warning,
    line_error,
    line_debug,
    line_verbosity,
    get_logger,
)

# XML I/O
from .xml_io import (
    NetworkXMLIO,
    NetworkXMLInfo,
    export_to_xml,
    import_from_xml,
    validate_xml,
    get_xml_info,
)

# JMT I/O
from .jmt_io import (
    JSIMGInfo,
    jmt2line,
    jsim2line,
    jmva2line,
    qn2jsimg,
)

# JMT Viewers
from .viewers import (
    get_jmt_path,
    jsimg_view,
    jsimw_view,
    jmva_view,
)

# Code generation
from .code_gen import (
    qn2python,
    line2python,
    sn2python,
    qn2java,
    lqn2java,
    line2java,
)

# Model converters
from .converters import (
    qn2line,
    line2qn,
    qn2lqn,
    lqn2qn,
    LQNTask,
    LQNProcessor,
    LQNModel,
    mapqn2renv,
    RandomEnvironmentModel,
    MMPP2Params,
)

# Model adapter
from .model_adapter import (
    ModelAdapter,
    DeaggregationInfo,
    TaggedModelResult,
    FESAggregationInfo,
    MMTResult,
    HTResult,
    tag_chain,
    aggregate_chains,
    remove_class,
    deaggregate_results,
    sort_forks,
    mmt,
    ht,
    paths,
    paths_cs,
    find_paths,
    find_paths_cs,
)

__all__ = [
    # Logging
    'VerboseLevel',
    'LineLogger',
    'line_printf',
    'line_warning',
    'line_error',
    'line_debug',
    'line_verbosity',
    'get_logger',
    # XML I/O
    'NetworkXMLIO',
    'NetworkXMLInfo',
    'export_to_xml',
    'import_from_xml',
    'validate_xml',
    'get_xml_info',
    # JMT I/O
    'JSIMGInfo',
    'jmt2line',
    'jsim2line',
    'jmva2line',
    'qn2jsimg',
    # JMT Viewers
    'get_jmt_path',
    'jsimg_view',
    'jsimw_view',
    'jmva_view',
    # Code generation
    'qn2python',
    'line2python',
    'sn2python',
    'qn2java',
    'lqn2java',
    'line2java',
    # Converters
    'qn2line',
    'line2qn',
    'qn2lqn',
    'lqn2qn',
    'LQNTask',
    'LQNProcessor',
    'LQNModel',
    'mapqn2renv',
    'RandomEnvironmentModel',
    'MMPP2Params',
    # Model adapter
    'ModelAdapter',
    'DeaggregationInfo',
    'TaggedModelResult',
    'FESAggregationInfo',
    'MMTResult',
    'HTResult',
    'tag_chain',
    'aggregate_chains',
    'remove_class',
    'deaggregate_results',
    'sort_forks',
    'mmt',
    'ht',
    'paths',
    'paths_cs',
    'find_paths',
    'find_paths_cs',
]
