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

#pragma once

/// LOCAL ADDITION (GEMS3K 2026-09-30): lets GEMS3K compile against an install with or without the two options below.
#define OPTIMA_LINESEARCH_STALL_ESCAPE 1
#define OPTIMA_LINESEARCH_REJECT_WORSE 1

namespace Optima {

/// The options for the line search minimization operation.
struct LineSearchOptions
{
    /// The tolerance in the minimization calculation during the line search operation.
    double tolerance = 1.0e-5;

    /// The maximum number of iterations during the minimization calculation in the line search operation.
    unsigned maxiterations = 20;

    /// The parameter that triggers line-search when current error is greater than initial error by a given factor (`Enew > factor*E0`).
    double trigger_when_current_error_is_greater_than_initial_error_by_factor = 1.0;

    /// The parameter that triggers line-search when current error is greater than previous error by a given factor (`Enew > factor*Eold`).
    double trigger_when_current_error_is_greater_than_previous_error_by_factor = 2.0;

    /// Whether ErrorControl may invoke the line search at all (LOCAL ADDITION,
    /// default false = upstream behaviour, where the call site is commented out).
    bool enabled = false;

    /// Whether the line search minimizes the UNMASKED residual norm instead of
    /// ResidualErrors::error() (LOCAL ADDITION, default false).
    /// error() zeroes the optimality residual of every unstable variable and of
    /// every basic variable sitting exactly on a bound, so it can be reduced by
    /// pushing variables onto their bounds rather than by making real progress -
    /// a perverse objective for a line search to minimize.
    bool use_unmasked_error = false;

    /// STALL ESCAPE (LOCAL ADDITION, GEMS3K 2026-09-30, default 0 = off): after this many CONSECUTIVE line-search executions
    /// that left the error unchanged (relative change <= stall_escape_tolerance), ErrorControl skips the line search once and
    /// keeps the full step. Measured on GEMS3K f_TestPNTDB: the strict trigger froze the error at exactly 2.29844 for 900
    /// iterations (AOP 1827 it); with 10 the solve takes 493 it, T-cement's line-search rescue is kept.
    std::size_t stall_escape_after = 0;

    /// Relative error change at or below which a line-search step counts as stalled. A looser value (1e-3, 1e-2) also caught
    /// GEMS3K Cu-Pourbaix's crawl but LOST T-cement's rescue at every setting tried.
    double stall_escape_tolerance = 1.0e-8;

    /// NON-MONOTONE WINDOW (LOCAL ADDITION, GEMS3K 2026-09-30, default 0 = off): compare the new error with the max of the
    /// last N pre-step errors instead of the previous one only. Helped GEMS3K Cu-Pourbaix (AOP 483 -> 110 it at 10) but was
    /// worse than the strict trigger on every corium phase diagram and lost T-cement - opt-in only.
    std::size_t nonmonotone_window = 0;

    /// REJECT IF WORSE (LOCAL ADDITION, GEMS3K 2026-09-30, default false = off): when a line search ends with an error NOT
    /// below the pre-step error, discard its result and keep the full step. Measured on GEMS3K (AOP, line search 1.5, stall
    /// escape 10): Cu-Pourbaix fired 341 line searches of which 329 ended ABOVE the pre-step error (a crawl the stall escape
    /// cannot see - each moves the error by ~0.3 %, in the wrong direction); with the rule 489 -> 149 it, j_Solvus 153 -> 104,
    /// same G; corium diagrams neutral on answers.
    bool reject_if_worse = false;
};

} // namespace Optima
