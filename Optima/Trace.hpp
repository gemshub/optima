// Optional solver trace, written when the REAKTORO_TRACE_FILE environment variable names a file.
#pragma once

#include <cstdio>
#include <vector>

namespace Optima {
namespace trace {

/// The trace file, or null when tracing is off. The variable is read only once.
FILE* file();

/// Counters shared between Optima and Reaktoro while tracing.
struct State
{
    long solve = 0;     ///< The number of solves started so far.
    int attempt = 1;    ///< 1 for the first attempt of a solve, 2 for its retry.
    long iter = 0;      ///< The iteration whose state is being evaluated.
    long fupdates = 0;  ///< The residual function updates in this solve.
    long fxevals = 0;   ///< The exact Hessian columns evaluated in this solve.
    std::vector<long> basis; ///< The basic variables last reported.
    bool auto_attempt = false; ///< When set, each Optima solve started within one traced solve counts as the next attempt.
    int solves_in_call = 0;    ///< The Optima solves started since beginSolve().
};

/// The shared counters.
State& state();

/// Start counting a new solve.
void beginSolve();

} // namespace trace
} // namespace Optima
