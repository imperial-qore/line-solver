/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_IO_QN2JSIMG_H
#define LINE_IO_QN2JSIMG_H

/**
 * @file
 * @ingroup line_io
 * Port of `matlab/src/io/QN2JSIMG.m`: write a network to a JMT `.jsimg` FILE
 * and return its path.
 *
 * `jmt_write_jsim` (jmt_writer.h) returns the DOCUMENT TEXT, because SolverJMT
 * stages it in a directory of its own. The reference's `QN2JSIMG` instead
 * returns `outputFileName`, the path written, and with no path given it writes
 * `<lineTempName('jsim')>/model.jsim`. This header supplies exactly that
 * contract on top of the text writer, so there is still one serializer.
 *
 * The temporary directory is NOT removed: the caller is handed its path to
 * open in JMT (`jsimgView`), as in the reference.
 */

#include <fstream>
#include <string>

#include "line/io/jmt_writer.h"
#include "line/lang/qn/network_builder.h"
#include "line/lang/qn/network_struct.h"
#include "line/util/error.h"
#include "line/util/tempdir.h"

namespace line {
namespace io {

/**
 * Port of `QN2JSIMG(model, outputFileName, options)` on a refreshed struct.
 *
 * @param sn    the refreshed struct (`model.getStruct(true)`)
 * @param path  the file to write; empty writes `model.jsim` in a fresh temp directory
 * @param opt   the simulation header controls, `JMTIO`'s properties
 * @return the path written
 */
template <class T>
std::string qn2jsimg(const qn::NetworkStruct<T>& sn, const std::string& path = std::string(),
                     const JmtWriteOptions& opt = JmtWriteOptions()) {
    const std::string out = path.empty() ? util::make_temp_dir("jsim") + "/model.jsim" : path;
    JmtWriteOptions o = opt;
    if (o.log_path.empty()) o.log_path = sn.log_path;
    const std::string text = jmt_write_jsim(sn, o);
    std::ofstream f(out.c_str(), std::ios::binary);
    if (!f) throw InputError("QN2JSIMG: cannot open '" + out + "' for writing");
    f.write(text.data(), static_cast<std::streamsize>(text.size()));
    f.close();
    if (!f) throw InputError("QN2JSIMG: failed to write '" + out + "'");
    return out;
}

/** Port of `QN2JSIMG(model, outputFileName, options)` on a model: refreshes it first. */
template <class T>
std::string qn2jsimg(qn::Network<T>& model, const std::string& path = std::string(),
                     const JmtWriteOptions& opt = JmtWriteOptions()) {
    return qn2jsimg(model.get_struct(), path, opt);
}

}  // namespace io
}  // namespace line

#endif  // LINE_IO_QN2JSIMG_H
