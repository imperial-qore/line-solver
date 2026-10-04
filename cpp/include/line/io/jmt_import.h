/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_JMT_IMPORT_H
#define LINE_IO_JMT_IMPORT_H

/**
 * @file
 * @ingroup line_io
 * The MATLAB-named JMT import entry points: `JSIM2LINE` and the `JMT2LINE`
 * dispatcher over `.jmva` / `.jsim` / `.jsimg` / `.jsimw`.
 *
 * `io::read_jsim` (jsim_reader.h) stays the implementation and keeps its name;
 * `jsim2line` adds only the reference's naming rule, `fileparts` of the
 * `<sim name=...>` attribute, so a model read from `model.jsimg` is called
 * `model` as it is in MATLAB and the JAR, where `read_jsim` keeps the raw
 * attribute.
 */

#include <iostream>
#include <string>

#include "line/io/jmva_reader.h"
#include "line/io/jsim_reader.h"
#include "line/lang/qn/network_builder.h"

namespace line {
namespace io {

/**
 * Port of `JSIM2LINE(filename, modelName)`.
 *
 * @param path  a `.jsim`, `.jsimg` or `.jsimw` file
 * @param name  the model name; empty takes `fileparts` of the document's `<sim name>`
 */
template <class T>
qn::Network<T> jsim2line(const std::string& path, const std::string& name = std::string()) {
    qn::Network<T> net = read_jsim<T>(path, name);
    if (name.empty()) net.raw_struct().name = jmva_detail::file_stem(net.raw_struct().name);
    return net;
}

namespace jmt_import_detail {

/** The extension without its dot, `fileparts`' third output minus the '.'. */
inline std::string extension(const std::string& path) {
    const std::size_t slash = path.find_last_of("/\\");
    const std::size_t dot = path.find_last_of('.');
    if (dot == std::string::npos || (slash != std::string::npos && dot < slash)) return "";
    return path.substr(dot + 1);
}

}  // namespace jmt_import_detail

/**
 * Port of `JMT2LINE(filename, modelName)`: dispatch on the extension.
 *
 * `.jmva` goes to `jmva2line`, `.jsim` / `.jsimg` / `.jsimw` to `jsim2line`.
 * Any other extension is read as JSIM after the reference's warning, which in
 * this port goes to stderr with the `[LINE] Warning:` prefix the io layer uses.
 * The match is case-sensitive, as `switch fext` is in MATLAB.
 */
template <class T>
qn::Network<T> jmt2line(const std::string& path, const std::string& name = std::string()) {
    const std::string ext = jmt_import_detail::extension(path);
    if (ext == "jmva") return jmva2line<T>(path, name);
    if (ext == "jsim" || ext == "jsimg" || ext == "jsimw") return jsim2line<T>(path, name);
    std::cerr << "[LINE] Warning: JMT2LINE: the file has unknown extension, trying to parse as a "
                 "JSIMG file."
              << std::endl;
    return jsim2line<T>(path, name);
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_JMT_IMPORT_H
