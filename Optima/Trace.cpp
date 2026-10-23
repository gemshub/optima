#include "Trace.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>

namespace Optima {
namespace trace {

FILE* file()
{
    static FILE* fp = []() -> FILE* {
        const char* fn = std::getenv("GEMS3K_OPTIMA_STEPTRACE_FILE");
        if(!fn) fn = std::getenv("REAKTORO_TRACE_FILE");   // the old name
        return fn ? std::fopen(fn, "a") : nullptr;
    }();
    return fp;
}

State& state()
{
    static State s;
    return s;
}

void beginSolve()
{
    auto& s = state();
    s.solve += 1;
    s.attempt = 1;
    s.iter = 0;
    s.fupdates = 0;
    s.fxevals = 0;
    s.basis.clear();
    s.solves_in_call = 0;
}

long countOutside(const double* x, const double* lo, const double* hi, long n)
{
    long c = 0;
    for(long i = 0; i < n; ++i)
        if(x[i] < lo[i] || x[i] > hi[i]) c++;
    return c;
}

void updateStepLimit(StepLimit& r, const double* x, const double* xo, const double* lo, const double* hi, long n, bool isP)
{
    for(long i = 0; i < n; ++i)
    {
        if(x[i] == xo[i]) continue;
        double b = 2.0, bd = 0.0;
        if(x[i] > hi[i] && xo[i] < hi[i]) { b = (hi[i] - xo[i])/(x[i] - xo[i]); bd = hi[i]; }
        else if(x[i] < lo[i] && xo[i] > lo[i]) { b = (lo[i] - xo[i])/(x[i] - xo[i]); bd = lo[i]; }
        if(b < r.beta) { r.beta = b; r.arg = isP ? -2 : i; r.xo = xo[i]; r.xnew = x[i]; r.bound = bd; }
    }
}

void writeNewton(const char* method, const double* dx, long nx, const double* dp, long np, const double* dw, long nw)
{
    auto* fp = file();
    if(!fp) return;
    long idu = 0;
    double mx = 0.0, mp = 0.0, mw = 0.0;
    for(long i = 0; i < nx; ++i)
        if(std::abs(dx[i]) > mx) { mx = std::abs(dx[i]); idu = i; }
    for(long i = 0; i < np; ++i) mp = std::max(mp, std::abs(dp[i]));
    for(long i = 0; i < nw; ++i) mw = std::max(mw, std::abs(dw[i]));
    fprintf(fp, "NEWTON solve=%ld attempt=%d it=%ld method=%s ndu_x=%.6e ndu_p=%.6e ndu_w=%.6e idu=%ld\n",
        state().solve, state().attempt, state().iter, method, mx, mp, mw, idu);
    fflush(fp);
}

void writeBeta(double betamin, const StepLimit& limit, long nclamped)
{
    auto* fp = file();
    if(!fp) return;
    fprintf(fp, "BETA solve=%ld attempt=%d it=%ld betamin=%.6e argmin=%ld xo=%.6e xnew=%.6e bound=%.6e n_clamped=%ld fix_and_accept=0\n",
        state().solve, state().attempt, state().iter, betamin, limit.arg, limit.xo, limit.xnew, limit.bound, nclamped);
    fflush(fp);
}

void writeBetaFixAndAccept(long nclamped)
{
    auto* fp = file();
    if(!fp) return;
    fprintf(fp, "BETA solve=%ld attempt=%d it=%ld betamin=1.000000e+00 argmin=-1 xo=0 xnew=0 bound=0 n_clamped=%ld fix_and_accept=1\n",
        state().solve, state().attempt, state().iter, nclamped);
    fflush(fp);
}

void Recorder::begin(const double* xlower, const double* xupper, long n)
{
    lo.assign(xlower, xlower + n);
    hi.assign(xupper, xupper + n);
    prevjb.clear();
    errors.clear();
}

void Recorder::iter(long it, double f, double err, double errx, double errp, double errw,
                    const std::ptrdiff_t* jb, long njb, const std::ptrdiff_t* jn, long njn, const double* x,
                    long nbs, long nl, long ns, long nlu, long nuu)
{
    auto* fp = file();
    if(!fp) return;
    long nbbound = 0, nnlow = 0, nnup = 0, jbchg = 0;
    for(long k = 0; k < njb; ++k) { const std::ptrdiff_t i = jb[k]; if(x[i] == lo[i] || x[i] == hi[i]) nbbound++; }
    for(long k = 0; k < njn; ++k) { const std::ptrdiff_t i = jn[k]; if(x[i] == lo[i]) nnlow++; else if(x[i] == hi[i]) nnup++; }
    std::vector<std::ptrdiff_t> cur(jb, jb + njb);
    for(auto i : cur) if(std::find(prevjb.begin(), prevjb.end(), i) == prevjb.end()) jbchg++;
    if(prevjb.empty()) jbchg = 0;
    prevjb = cur;
    errors.push_back(err);
    const auto& ts = state();
    fprintf(fp, "ITER solve=%ld attempt=%d it=%ld f=%.12e err=%.6e errx=%.6e errp=%.6e errw=%.6e nb=%ld nbs=%ld nl=%ld ns=%ld nlu=%ld nuu=%ld nn_lower=%ld nn_upper=%ld nb_bound=%ld jb_changed=%ld F_updates=%ld fx_evals=%ld\n",
        ts.solve, ts.attempt, it, f, err, errx, errp, errw, njb, nbs, nl, ns, nlu, nuu, nnlow, nnup, nbbound, jbchg, ts.fupdates, ts.fxevals);
    fflush(fp);
}

void Recorder::stop(long it, const char* reason, double err, double tol, unsigned maxiters, bool check)
{
    auto* fp = file();
    if(!fp) return;
    long i6 = -1, i8 = -1;
    for(size_t i = 0; i < errors.size(); ++i)
    {
        if(i6 < 0 && errors[i] < 1e-6) i6 = (long)i;
        if(i8 < 0 && errors[i] < 1e-8) i8 = (long)i;
    }
    fprintf(fp, "STOP solve=%ld attempt=%d it=%ld reason=%s err=%.6e tol_used=%.3e maxiters_used=%u check_used=%d first_it_below_1e-6=%ld first_it_below_1e-8=%ld\n",
        state().solve, state().attempt, it, reason, err, tol, maxiters, check ? 1 : 0, i6, i8);
    fflush(fp);
}

} // namespace trace
} // namespace Optima
