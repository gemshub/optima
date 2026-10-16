#include "Trace.hpp"

#include <cstdlib>

namespace Optima {
namespace trace {

FILE* file()
{
    static FILE* fp = []() -> FILE* {
        const char* fn = std::getenv("REAKTORO_TRACE_FILE");
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

} // namespace trace
} // namespace Optima
