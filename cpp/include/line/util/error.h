/*
 * Copyright (c) 2012-2026, QORE Lab, Imperial College London
 * All rights reserved.
 */
#ifndef LINE_UTIL_ERROR_H
#define LINE_UTIL_ERROR_H

#include <stdexcept>
#include <string>

namespace line {

/**
 * Base error for the multiprecision C++ port. Nothing in the library calls
 * exit() or prints to stderr on failure: the CLI catches these at main and the
 * Python binding maps them to exceptions, so the same code path serves both.
 */
class Error : public std::runtime_error {
public:
    explicit Error(const std::string& what) : std::runtime_error(what) {}
};

/** Malformed or inconsistent input (dimensions, negative populations, ...). */
class InputError : public Error {
public:
    explicit InputError(const std::string& what) : Error(what) {}
};

/** The algorithm cannot proceed on this instance (singular matrix, ...). */
class NumericError : public Error {
public:
    explicit NumericError(const std::string& what) : Error(what) {}
};

/** Requested feature or arithmetic mode is not ported yet. */
class UnsupportedError : public Error {
public:
    explicit UnsupportedError(const std::string& what) : Error(what) {}
};

}  // namespace line

#endif  // LINE_UTIL_ERROR_H
