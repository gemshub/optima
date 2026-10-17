#include <cmath>
// Optima is a C++ library for solving linear and non-linear constrained optimization problems.
//
// Copyright © 2020-2024 Allan Leal
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.

#include <cstdio>
#include <cstdlib>
#include "BacktrackSearch.hpp"

// Optima includes
#include <Optima/Exception.hpp>
#include <Optima/Trace.hpp>

namespace Optima {

using std::min;
using std::max;
using std::abs;
using std::greater;

struct BacktrackSearch::Impl
{
    BacktrackSearchOptions options;  ///< The options for the backtrack search operation.
    MasterDims dims;                 ///< The dimensions of the master variables.
    MasterVector unew;               ///< The state of u = (x, p, y, z) right-after Newton step without any correction.
    Vector xlower;                   ///< The lower bounds for x.
    Vector xupper;                   ///< The upper bounds for x.
    Vector plower;                   ///< The lower bounds for p.
    Vector pupper;                   ///< The upper bounds for p.
    Vector betas;                    ///< The beta factors for x and p

    Impl()
    {
    }

    auto initialize(const MasterProblem& problem) -> void
    {
        dims   = problem.dims;
        xlower = problem.xlower;
        xupper = problem.xupper;
        plower = problem.plower;
        pupper = problem.pupper;
        unew.resize(dims);
    }

    auto execute(MasterVectorView uo, MasterVectorRef u, ResidualFunction& F, ResidualErrors& E) -> void
    {
        auto const& xo = uo.x;
        auto const& po = uo.p;
        auto const& x = u.x;
        auto const& p = u.p;
        assert((xupper.array() >= xo.array()).all());
        assert((xlower.array() <= xo.array()).all());
        assert((pupper.array() >= po.array()).all());
        assert((plower.array() <= po.array()).all());

        if(options.apply_min_max_fix_and_accept)
        {
            const long nclamp = trace::file() ? trace::countOutside(u.x.data(), xlower.data(), xupper.data(), (long)dims.nx) : 0;
            u.x.noalias() = min(max(u.x, xlower), xupper);
            u.p.noalias() = min(max(u.p, plower), pupper);
            trace::writeBetaFixAndAccept(nclamp);
            return;
        }

        auto betamin = 1.0;

        int probeArg = -1; double probeXo = 0., probeXn = 0., probeBound = 0.;

        for(auto i = 0; i < dims.nx; ++i)
        {
            if(x[i] == xo[i])
                continue;
            if(x[i] > xupper[i] && xo[i] < xupper[i]) {
                const auto b = (xupper[i] - xo[i])/(x[i] - xo[i]);
                if(b < betamin) { betamin = b; probeArg = i; probeXo = xo[i]; probeXn = x[i]; probeBound = xupper[i]; }
            }
            else if(x[i] < xlower[i] && xo[i] > xlower[i]) {
                const auto b = (xlower[i] - xo[i])/(x[i] - xo[i]);
                if(b < betamin) { betamin = b; probeArg = i; probeXo = xo[i]; probeXn = x[i]; probeBound = xlower[i]; }
            }
        }

        for(auto i = 0; i < dims.np; ++i)
        {
            if(p[i] == po[i])
                continue;
            if(p[i] > pupper[i] && po[i] < pupper[i])
                { const auto b = (pupper[i] - po[i])/(p[i] - po[i]); if(b < betamin) { betamin = b; probeArg = -2; probeXo = po[i]; probeXn = p[i]; probeBound = pupper[i]; } }
            else if(p[i] < plower[i] && po[i] > plower[i])
                { const auto b = (plower[i] - po[i])/(p[i] - po[i]); if(b < betamin) { betamin = b; probeArg = -2; probeXo = po[i]; probeXn = p[i]; probeBound = plower[i]; } }
        }

        // Diagnostic: if the OPTIMA_BETA_PROBE environment variable names a file, write
        // the final betamin and the variable that sets it there on every call
        // (argmin -1: no limit, -2: a p variable, xo/xnew/bound then refer to it).
        {
            static FILE* bpf = [] {
                const char* bp = std::getenv("OPTIMA_BETA_PROBE");
                return bp ? fopen(bp, "w") : nullptr;
            }();
            if(bpf) {
                fprintf(bpf, "betamin %.6e argmin %d xo %.6e xnew %.6e bound %.6e\n",
                        betamin, probeArg, probeXo, probeXn, probeBound);
            }
        }

        // Find which variable sets betamin (same rule as above), only when tracing.
        trace::StepLimit tlimit;
        if(trace::file())
        {
            trace::updateStepLimit(tlimit, x.data(), xo.data(), xlower.data(), xupper.data(), (long)dims.nx, false);
            trace::updateStepLimit(tlimit, p.data(), po.data(), plower.data(), pupper.data(), (long)dims.np, true);
        }

        // Keep a variable whose step is zero exactly where it was: the blend below can
        // round it one unit away, which frees a variable held on its bound.
        const Vector xstep = u.x - uo.x;
        const Vector pstep = u.p - uo.p;

        u = uo*(1 - betamin) + betamin*u;

        for(auto i = 0; i < dims.nx; ++i) if(xstep[i] == 0.0) u.x[i] = uo.x[i];
        for(auto i = 0; i < dims.np; ++i) if(pstep[i] == 0.0) u.p[i] = uo.p[i];

        // A variable whose full step crossed a bound and that the blend left within a few units
        // of it is put exactly on the bound: one unit short, it is never held there and every
        // later step is near zero.
        auto snap = [](double v, double target, double lo, double hi) {
            auto near = [](double a, double bd) { return std::abs(a - bd) <= 4.0*std::abs(std::nextafter(bd, 2.0*bd + 1.0) - bd); };
            if(target < lo && near(v, lo)) return lo;
            if(target > hi && near(v, hi)) return hi;
            return v; };
        for(auto i = 0; i < dims.nx; ++i) u.x[i] = snap(u.x[i], uo.x[i] + xstep[i], xlower[i], xupper[i]);
        for(auto i = 0; i < dims.np; ++i) u.p[i] = snap(u.p[i], uo.p[i] + pstep[i], plower[i], pupper[i]);

        const long nclamp = trace::file() ? trace::countOutside(u.x.data(), xlower.data(), xupper.data(), (long)dims.nx) : 0;

        u.x.noalias() = min(max(u.x, xlower), xupper);
        u.p.noalias() = min(max(u.p, plower), pupper);

        trace::writeBeta(betamin, tlimit, nclamp);
    }
};

BacktrackSearch::BacktrackSearch()
: pimpl(new Impl())
{}

BacktrackSearch::BacktrackSearch(const BacktrackSearch& other)
: pimpl(new Impl(*other.pimpl))
{}

BacktrackSearch::~BacktrackSearch()
{}

auto BacktrackSearch::operator=(BacktrackSearch other) -> BacktrackSearch&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto BacktrackSearch::setOptions(const BacktrackSearchOptions& options) -> void
{
    pimpl->options = options;
}

auto BacktrackSearch::initialize(const MasterProblem& problem) -> void
{
    pimpl->initialize(problem);
}

auto BacktrackSearch::execute(MasterVectorView uo, MasterVectorRef u, ResidualFunction& F, ResidualErrors& E) -> void
{
    pimpl->execute(uo, u, F, E);
}

} // namespace Optima
