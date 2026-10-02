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
            u.x.noalias() = min(max(u.x, xlower), xupper);
            u.p.noalias() = min(max(u.p, plower), pupper);
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

        u = uo*(1 - betamin) + betamin*u;

        u.x.noalias() = min(max(u.x, xlower), xupper);
        u.p.noalias() = min(max(u.p, plower), pupper);
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
