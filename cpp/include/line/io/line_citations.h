// Copyright (c) 2012-2026, QORE Lab, Imperial College London
// All rights reserved.
#ifndef LINE_IO_LINE_CITATIONS_H
#define LINE_IO_LINE_CITATIONS_H

/**
 * @file
 * Bibliographic references for the algorithms a run used.
 *
 * Port of MATLAB `matlab/src/io/line_citations.m`, the JAR
 * `jline.io.LineCitations` and python `line_solver/api/io/citations.py`. The
 * registry is the manual's method-to-citation table (`doc/latex/manual.tex`)
 * and its bibliography; keep the two in step.
 *
 * ATTRIBUTION IN LINE IS PULL-BASED. Nothing is printed during a solve; a user
 * asks for the references when writing the run up. A method that reaches a user
 * without an entry here is a silent loss of attribution, so a new or renamed
 * solver method, analyzer, transformation or approximation adds its entry in
 * the SAME change, in all four codebases.
 *
 * METHOD NAME LOOKUP, identical to the other three: a method name is lowercased and
 * trimmed; an exact match wins; otherwise a family-qualified method name falls back
 * to the bare method name after its first dot, so `nc.mva` finds `mva` when
 * `nc.mva` is absent. Unknown method names are ignored, so a caller may pass whatever
 * it knows about a run, including the `default/<actual>` compound a dispatching
 * solver reports (split it on '/' and pass both halves). Results are
 * de-duplicated BY BIBLIOGRAPHY KEY, since one paper is often reachable through
 * several method names, and keep the order of first appearance.
 */

#include <map>
#include <string>
#include <vector>

namespace line {
namespace io {

/** One bibliography entry. */
struct Citation {
    /** Bibliography key, as used in doc/latex/biblio.bib. */
    std::string key;
    /** Short reference: author, title, venue, year. */
    std::string ref;
    /** One line saying which part of the solution process it covers. */
    std::string covers;
};

/** The whole method name -> reference table, built once. */
const std::map<std::string, Citation>& citation_registry();

/** The references for METHOD NAMES, de-duplicated by key, in order of appearance. */
std::vector<Citation> line_citations(const std::vector<std::string>& method_names);

/** Convenience for the `default/<actual>` compound a solver reports. */
std::vector<Citation> line_citations_for_method(const std::string& method);

}  // namespace io
}  // namespace line

#endif  // LINE_IO_LINE_CITATIONS_H
