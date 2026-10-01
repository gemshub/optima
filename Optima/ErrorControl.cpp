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

#include "ErrorControl.hpp"
#include <vector>
#include <algorithm>
#include <cmath>

// Optima includes
#include <Optima/Exception.hpp>
#include <Optima/BacktrackSearch.hpp>
#include <Optima/ErrorStatus.hpp>
#include <Optima/LineSearch.hpp>

namespace Optima {

struct ErrorControl::Impl
{
    ErrorStatus errorstatus;         ///< The current error status of the calculation.
    BacktrackSearch backtracksearch; ///< The backtrack algorithm to correct steps producing infinity errors.
    LineSearch linesearch;           ///< The line-search algorithm to correct steps producing significant large errors.

    Impl()
    {}

    ErrorControlOptions options; ///< LOCAL ADDITION: kept so execute() can read the line-search trigger.

    std::vector<double> errhist;     ///< LOCAL ADDITION: pre-step errors for linesearch.nonmonotone_window
    std::size_t stallCount = 0;      ///< LOCAL ADDITION: consecutive zero-progress line searches (linesearch.stall_escape_after)

    auto setOptions(const ErrorControlOptions& opts) -> void
    {
        options = opts;
        errorstatus.setOptions(options.errorstatus);
        backtracksearch.setOptions(options.backtracksearch);
        linesearch.setOptions(options.linesearch);
    }

    auto initialize(const MasterProblem& problem) -> void
    {
        errhist.clear();
        stallCount = 0;
        errorstatus.initialize();
        backtracksearch.initialize(problem);
        linesearch.initialize(problem);
    }

    auto execute(MasterVectorView uo, MasterVectorRef u, ResidualFunction& F, ResidualErrors& E) -> void
    {
        // LOCAL ADDITION: when the line search is asked to minimize the
        // unmasked residual norm, the TRIGGER must read the same quantity.
        // Measured (GEMS3K/CLAUDE.md 2026-08-25): with the trigger on
        // E.error(), no project in the GEMS3K suite ever fires the line search
        // at any factor >= 1, because the masked error essentially never grows
        // between iterations - it can be reduced by pushing variables onto
        // their bounds, which is exactly what makes it a poor merit function.
        const auto use_raw = options.linesearch.enabled && options.linesearch.use_unmasked_error;

        const auto error_prev = use_raw ? E.errorRaw() : E.error();

        backtracksearch.execute(uo, u, F, E);

        // LOCAL FIX (GEMS3K plan v5 s139.6, 2026-09-28, owner-approved): nothing between error_prev and here
        // updates E - BacktrackSearch changes only u, and MasterSolver::step() calls F.update/E.update only
        // AFTER this function - so error_new used to equal error_prev and the trigger compared a number with
        // itself (never fires at a factor >= 1; fires every step with a bare >=). When the line search is
        // enabled, evaluate the error at the new point so the comparison is real. Disabled (the default),
        // nothing here runs and the solver is unchanged. The extra evaluation is not side-effect free in
        // GEMS3K's objective, which is why it is confined to the enabled case.
        if(options.linesearch.enabled) { F.update(u); E.update(u, F); }

        const auto error_new = use_raw ? E.errorRaw() : E.error();

        // LOCAL ADDITION. Upstream ships this call commented out; LineSearch.cpp
        // is nonetheless fully implemented, and LineSearchOptions carries two
        // trigger-factor fields that nothing reads - which suggests the intended
        // trigger was the factor below, not the bare `error_new >= error_prev`
        // the commented-out line uses. Enabling that bare form was measured to be
        // a severe regression (GEMS3K/CLAUDE.md 2026-08-25), so this stays opt-in
        // (options.linesearch.enabled, default false = upstream behaviour) and
        // uses the factor.
        // LOCAL ADDITION (GEMS3K 2026-09-30): optional non-monotone reference - max of the last N pre-step errors.
        double error_ref = error_prev;
        const auto nmN = options.linesearch.nonmonotone_window;
        if( nmN > 0 )
        {
            errhist.push_back(error_prev);
            if( errhist.size() > nmN ) errhist.erase(errhist.begin());
            for( double e : errhist ) error_ref = std::max(error_ref, e);
        }
        const auto stallK = options.linesearch.stall_escape_after;
        if( options.linesearch.enabled &&
            error_new > options.linesearch.trigger_when_current_error_is_greater_than_previous_error_by_factor * error_ref )
        {
            // LOCAL ADDITION (GEMS3K 2026-09-30): stall escape - after stallK consecutive zero-progress line searches keep
            // the full step once (F and E are already evaluated at u above).
            if( stallK > 0 && stallCount >= stallK )
            {
                stallCount = 0;
                return;
            }
            // LOCAL ADDITION (GEMS3K 2026-09-30): reject_if_worse - a line search that ends at or above the pre-step error
            // is discarded in favour of the full step.
            MasterVector ufull;
            if( options.linesearch.reject_if_worse ) ufull = u;
            linesearch.execute(uo, u, F, E);
            if( options.linesearch.reject_if_worse && ( use_raw ? E.errorRaw() : E.error() ) >= error_prev )
            {
                u = ufull;
                F.update(u);
                E.update(u, F);
            }
            if( stallK > 0 )
            {
                const auto error_ls = use_raw ? E.errorRaw() : E.error();
                const double rel = std::abs(error_ls - error_prev) / std::max(std::abs(error_prev), 1e-300);
                stallCount = ( rel <= options.linesearch.stall_escape_tolerance ) ? stallCount + 1 : 0;
            }
        }
        else if( stallK > 0 ) stallCount = 0;
    }
};

ErrorControl::ErrorControl()
: pimpl(new Impl())
{}

ErrorControl::ErrorControl(const ErrorControl& other)
: pimpl(new Impl(*other.pimpl))
{}

ErrorControl::~ErrorControl()
{}

auto ErrorControl::operator=(ErrorControl other) -> ErrorControl&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto ErrorControl::setOptions(const ErrorControlOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto ErrorControl::initialize(const MasterProblem& problem) -> void
{
    pimpl->initialize(problem);
}

auto ErrorControl::execute(MasterVectorView uo, MasterVectorRef u, ResidualFunction& F, ResidualErrors& E) -> void
{
    pimpl->execute(uo, u, F, E);
}

} // namespace Optima
