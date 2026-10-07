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

#include <cstddef>

/// Defined when LineSearchOptions has stall_escape_after and reject_if_worse, so client code can test for them.
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

    /// Whether ErrorControl runs the line search (default false = never).
    /// In plain words: allow the solver to try a shorter step when a full step makes things worse.
    bool enabled = false;

    /// Whether the line search minimizes ResidualErrors::errorRaw() instead of error() (default false).
    /// In plain words: judge progress by the full error, including the parts error() hides
    /// for variables sitting on a bound.
    bool use_unmasked_error = false;

    /// After this many consecutive line searches that leave the error unchanged
    /// (relative change <= stall_escape_tolerance), skip the line search once and keep the
    /// full step (default 0 = off).
    /// In plain words: if the shorter steps keep getting nowhere, take one full step to break out.
    std::size_t stall_escape_after = 0;

    /// Relative error change at or below which a line search counts as making no progress.
    double stall_escape_tolerance = 1.0e-8;

    /// If a line search ends with an error not below the pre-step error, discard it and keep
    /// the full step (default false = off).
    /// In plain words: if the shorter step did not help, use the full step after all.
    bool reject_if_worse = false;
};

} // namespace Optima
