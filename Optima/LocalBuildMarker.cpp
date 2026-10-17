// Prints a line to stderr when the library is loaded, to show that this
// modified Optima build is the one in use and not a stock one.
// Only in builds made with -DOPTIMA_LOAD_MARKER=ON; release builds stay quiet.
#ifdef OPTIMA_LOAD_MARKER
#include <cstdio>

namespace {
struct LocalBuildMarker
{
    LocalBuildMarker()
    {
        std::fprintf( stderr,
            "[Optima] LOCAL MODIFIED BUILD loaded (gemshub fork) - if you expected stock Optima, "
            "check your RPATH/CMAKE_PREFIX_PATH ordering.\n" );
    }
};
LocalBuildMarker g_localBuildMarker;
}
#endif
