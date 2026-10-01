// Prints a line to stderr when the library is loaded, to show that this
// modified Optima build is the one in use and not a stock one.
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
