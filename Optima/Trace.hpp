// Optional solver trace, written when the GEMS3K_OPTIMA_STEPTRACE_FILE environment variable (or the old
// name REAKTORO_TRACE_FILE) names a file.
#pragma once

#include <cstddef>
#include <cstdio>
#include <vector>

namespace Optima {
namespace trace {

/// The trace file, or null when tracing is off. The variable is read only once.
FILE* file();

/// Counters shared between Optima and the program that calls it, while tracing.
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
    long snaps = 0;            ///< Variables put on their lower bound by snap_leaving_after (counted always).
    long lupartial = 0;        ///< Solves done with the partial-pivot LU (counted always).
    long lufallbacks = 0;      ///< Partial-pivot solves redone with full pivoting (counted always).
};

/// The shared counters.
State& state();

/// Start counting a new solve.
void beginSolve();

/// The number of entries of x outside [lo, hi].
long countOutside(const double* x, const double* lo, const double* hi, long n);

/// The variable that limits a step (arg = -1 when none does, -2 for a p variable) and the step length it allows.
struct StepLimit { long arg = -1; double xo = 0, xnew = 0, bound = 0, beta = 1.0; };

/// Update the step limit with the variables of one block (x, or p when isP) stepping from xo to x inside [lo, hi].
void updateStepLimit(StepLimit& r, const double* x, const double* xo, const double* lo, const double* hi, long n, bool isP);

/// Write a NEWTON record: the method and the largest changes of x (with its index), p and w.
void writeNewton(const char* method, const double* dx, long nx, const double* dp, long np, const double* dw, long nw);

/// Write a BETA record: the step length and what limited it.
void writeBeta(double betamin, const StepLimit& limit, long nclamped);

/// Write the BETA record of a full step clipped to the bounds.
void writeBetaFixAndAccept(long nclamped);

/// The ITER and STOP records of one solve.
class Recorder
{
public:
    /// Start a solve with these bounds of x.
    void begin(const double* xlower, const double* xupper, long n);

    /// Write one ITER record.
    void iter(long it, double f, double err, double errx, double errp, double errw,
              const std::ptrdiff_t* jb, long njb, const std::ptrdiff_t* jn, long njn, const double* x,
              long nbs, long nl, long ns, long nlu, long nuu);

    /// Write the STOP record: why the loop ended and when the error first passed 1e-6 and 1e-8.
    void stop(long it, const char* reason, double err, double tol, unsigned maxiters, bool check);

private:
    std::vector<double> lo, hi;     ///< The bounds of x.
    std::vector<std::ptrdiff_t> prevjb;       ///< The basic variables at the previous iteration.
    std::vector<double> errors;     ///< The error at each iteration.
};

} // namespace trace
} // namespace Optima
