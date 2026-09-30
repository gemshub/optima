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
    Vector betasx;                   ///< Per-variable beta factors for x (ORCHESTRA mode).
    Vector betasp;                   ///< Per-variable beta factors for p (ORCHESTRA mode).

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

    /// Per-variable fraction-to-the-boundary factor for a single component.
    /// Returns 1.0 when the component does not leave its box (and, when a
    /// relative trust region is active, does not leave that either).
    auto componentBeta(double xo_i, double x_i, double lo_i, double hi_i, bool apply_ratio) const -> double
    {
        if(x_i == xo_i)
            return 1.0;
        auto beta = 1.0;
        if(x_i > hi_i && xo_i < hi_i)
            beta = min(beta, (hi_i - xo_i)/(x_i - xo_i));
        else if(x_i < lo_i && xo_i > lo_i)
            beta = min(beta, (lo_i - xo_i)/(x_i - xo_i));
        if(apply_ratio && xo_i > 0.0)
        {
            const auto r = options.max_step_ratio;
            const auto rlo = xo_i / r;
            const auto rhi = xo_i * r;
            if(x_i < rlo)
                beta = min(beta, (rlo - xo_i)/(x_i - xo_i));
            else if(x_i > rhi)
                beta = min(beta, (rhi - xo_i)/(x_i - xo_i));
        }
        return beta;
    }

    /// ORCHESTRA-style step control - see BacktrackSearchOptions::min_common_beta.
    /// Two parts, both load-bearing: (1) a variable whose own beta is at or
    /// below the threshold does not set the shared betamin, and (2) it is still
    /// clamped to its own beta, i.e. every variable moves by
    /// min(own_beta, betamin).
    auto executeOrchestra(MasterVectorView uo, MasterVectorRef u) -> void
    {
        auto const& xo = uo.x;
        auto const& po = uo.p;
        auto const& wo = uo.w;

        const auto fmin = options.min_common_beta;
        const auto apply_ratio = options.max_step_ratio > 1.0;

        betasx.resize(dims.nx);
        betasp.resize(dims.np);

        auto betamin = 1.0;

        for(auto i = 0; i < dims.nx; ++i)
        {
            const auto b = componentBeta(xo[i], u.x[i], xlower[i], xupper[i], apply_ratio);
            betasx[i] = b;
            if(b > fmin)
                betamin = min(betamin, b);
        }

        for(auto i = 0; i < dims.np; ++i)
        {
            const auto b = componentBeta(po[i], u.p[i], plower[i], pupper[i], false);
            betasp[i] = b;
            if(b > fmin)
                betamin = min(betamin, b);
        }

        for(auto i = 0; i < dims.nx; ++i)
        {
            const auto s = min(betasx[i], betamin);
            u.x[i] = xo[i]*(1 - s) + s*u.x[i];
        }

        for(auto i = 0; i < dims.np; ++i)
        {
            const auto s = min(betasp[i], betamin);
            u.p[i] = po[i]*(1 - s) + s*u.p[i];
        }

        u.w.noalias() = wo*(1 - betamin) + betamin*u.w;

        u.x.noalias() = min(max(u.x, xlower), xupper);
        u.p.noalias() = min(max(u.p, plower), pupper);
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

        if(options.min_common_beta > 0.0)
        {
            executeOrchestra(uo, u);
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

        // Diagnostic: which variable sets the shared `betamin`, and how small it gets.
        // Off unless OPTIMA_BETA_PROBE names a file; the env lookup and fopen happen once,
        // so the disabled path is a single null-pointer test per call.
        //
        // Measured on GEMS3K's seawater case (2026-09-01), 10 000 iterations: betamin NEVER
        // reaches 1.0, median 1.7e-16, and 97 % of the sub-1e-6 values are set by a variable
        // already sitting AT its lower bound whose Newton step points negative. One such
        // variable annihilates the step for all of them, because betamin is a shared scalar.
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

        for(auto i = 0; i < dims.np; ++i)
        {
            if(p[i] == po[i])
                continue;
            if(p[i] > pupper[i] && po[i] < pupper[i])
                betamin = min(betamin, (pupper[i] - po[i])/(p[i] - po[i]));
            else if(p[i] < plower[i] && po[i] > plower[i])
                betamin = min(betamin, (plower[i] - po[i])/(p[i] - po[i]));
        }

        if(options.max_step_ratio > 1.0)
        {
            const auto r = options.max_step_ratio;
            for(auto i = 0; i < dims.nx; ++i)
            {
                if(x[i] == xo[i] || xo[i] <= 0.0)
                    continue;
                const auto lo = xo[i] / r;
                const auto hi = xo[i] * r;
                if(x[i] < lo)
                    betamin = min(betamin, (lo - xo[i])/(x[i] - xo[i]));
                else if(x[i] > hi)
                    betamin = min(betamin, (hi - xo[i])/(x[i] - xo[i]));
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
