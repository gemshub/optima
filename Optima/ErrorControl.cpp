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

    ErrorControlOptions options; ///< The options for the error control.

    std::size_t stallCount = 0;      ///< Consecutive line searches that did not change the error (linesearch.stall_escape_after).

    auto setOptions(const ErrorControlOptions& opts) -> void
    {
        options = opts;
        errorstatus.setOptions(options.errorstatus);
        backtracksearch.setOptions(options.backtracksearch);
        linesearch.setOptions(options.linesearch);
    }

    auto initialize(const MasterProblem& problem) -> void
    {
        stallCount = 0;
        errorstatus.initialize();
        backtracksearch.initialize(problem);
        linesearch.initialize(problem);
    }

    auto execute(MasterVectorView uo, MasterVectorRef u, ResidualFunction& F, ResidualErrors& E) -> void
    {
        // The line-search trigger uses the same error measure the line search minimizes.
        const auto use_raw = options.linesearch.enabled && options.linesearch.use_unmasked_error;

        const auto error_prev = use_raw ? E.errorRaw() : E.error();

        backtracksearch.execute(uo, u, F, E);

        // With the line search enabled, evaluate the error at the new point; otherwise
        // error_new would still be the pre-step error.
        if(options.linesearch.enabled) { F.update(u); E.update(u, F); }

        const auto error_new = use_raw ? E.errorRaw() : E.error();

        // Run the line search when the new error exceeds the pre-step error by the trigger factor.
        // In plain words: if the full step made the error clearly worse, try a shorter step.
        const auto stallK = options.linesearch.stall_escape_after;
        if( options.linesearch.enabled &&
            error_new > options.linesearch.trigger_when_current_error_is_greater_than_previous_error_by_factor * error_prev )
        {
            // After stallK consecutive line searches that did not change the error, keep the full step once.
            if( stallK > 0 && stallCount >= stallK )
            {
                stallCount = 0;
                return;
            }
            // reject_if_worse: if the line search ends at or above the pre-step error, keep the full step instead.
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
