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

namespace Optima {

/// The options for the backtrack search operations.
struct BacktrackSearchOptions
{
    /// The flag that indicates if a simpler approach for fixing out-of-bounds
    /// variables should be used. Setting this option to true will cause the
    /// updated variables `x` and `p` after a Newton step to be fixed with `x =
    /// min(xlower, max(x, xupper))` and `p = min(plower, max(p, pupper))` and
    /// immediately accepted. If line-search is to be performed afterwards, it
    /// is possible that the new direction vector from previous to updated `u`
    /// vector is not a descent direction, since the min-max fix alters its
    /// orientation.
    bool apply_min_max_fix_and_accept = false;

    /// Relative step limit for `x`. When `r > 1`, the whole step is scaled down
    /// so that no `x[i]` with `xo[i] > 0` leaves `[xo[i]/r, xo[i]*r]`; the step
    /// direction is kept. 0.0 (default) = off.
    /// In plain words: no amount may grow or shrink by more than a factor r in one step.
    double max_step_ratio = 0.0;

    /// Per-variable step control. When `f > 0`, a variable whose own
    /// fraction-to-the-boundary factor is <= `f` does not limit the shared
    /// factor `betamin`; each variable then moves by `min(own_beta, betamin)`,
    /// and the duals `w` by `betamin`. The step direction is not kept, so the
    /// linear equality constraints may be violated after the step.
    /// 0.0 (default) = off.
    /// In plain words: a variable stuck at its limit no longer shortens the step of all the others.
    double min_common_beta = 0.0;
};

} // namespace Optima
