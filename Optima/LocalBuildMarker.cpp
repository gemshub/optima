// Not part of upstream Optima. Added 2026-08-24 while integrating this
// modified checkout into GEMS3K (see BacktrackSearchOptions::max_step_ratio
// and MasterSolver.cpp's convergence-check reordering).
//
// Purpose: this checkout has both a stock, conda-packaged Optima (used by
// Reaktoro's own prebuilt Python package, and previously cached in a build
// tree's CMakeCache.txt / an ad-hoc test binary's own RPATH) and this
// modified one coexisting on the same machine, resolvable from more than
// one CMAKE_PREFIX_PATH entry. A real, previously-hit failure mode (see
// GEMS3K's CLAUDE.md, "Testing/debugging gotcha") is a consumer silently
// linking/loading the WRONG one - which, for a change like
// max_step_ratio that only exists here, doesn't fail to compile or link
// (the GEMS3K side compiles fine against this checkout's headers) - it
// produces an ABI mismatch at runtime instead (Options's layout differs
// between the two builds), silently reading/writing past the true struct
// bounds. That is far worse than a clean failure, and exactly what
// happened once already this session before this marker existed - the
// first "it works" result on Resources/gems3k/j_Flowline_G_series1_...
// was undefined behavior from the wrong library being loaded, not a real
// success.
//
// This unconditional, load-time (static-init, so it fires the moment the
// shared library is mapped into a process - no API call required) stderr
// print is the cheapest possible unambiguous confirmation that a given
// process really did load THIS build, not a stock one - check for this
// exact line in any tool's stderr output before trusting a result that
// depends on this checkout's own additions.
#include <cstdio>

namespace {
struct LocalBuildMarker
{
    LocalBuildMarker()
    {
        std::fprintf( stderr,
            "[Optima] LOCAL MODIFIED BUILD loaded (checkout: /home/dmiron/git/hub/optima, "
            "carries BacktrackSearchOptions::max_step_ratio + MasterSolver convergence-order fix, "
            "2026-08-24) - if you expected stock Optima, check your RPATH/CMAKE_PREFIX_PATH ordering.\n" );
    }
};
LocalBuildMarker g_localBuildMarker;
}
