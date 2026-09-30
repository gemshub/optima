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

    /// Optional relative trust-region cap on how much any component of `x`
    /// may grow or shrink in a single (uncorrected) Newton step, applied via
    /// the same shared betamin scale-down already used for box-bound
    /// clipping below (i.e. it scales the WHOLE step, preserving the Newton
    /// direction, rather than clamping each variable independently).
    /// Disabled (0.0, the default) reproduces the exact prior behavior with
    /// zero effect on any existing consumer. When set to a value `r > 1`,
    /// no `x[i]` with `xo[i] > 0` is allowed to move outside
    /// `[xo[i]/r, xo[i]*r]` in one step. Added to guard against a strictly-
    /// positive interior variable (no active box bound nearby) being driven
    /// to its numerical floor by an undamped Newton step in one iteration -
    /// box-bound backtracking alone does not catch this, since the bound
    /// itself may be many orders of magnitude below the variable's current
    /// value.
    double max_step_ratio = 0.0;

    /// ORCHESTRA-style outlier-decoupled step control (default 0.0 = disabled,
    /// exactly reproducing the prior shared-betamin behavior).
    ///
    /// Ported from ORCHESTRA's `UnEqGroup::adaptEstimations()`
    /// (orchestra_cpp/UnEqGroup.cpp:564-608, `minimumfactor = 1e-5`). The
    /// shared backtracking factor `betamin` is normally the minimum over every
    /// variable's own fraction-to-the-boundary factor, so a single variable
    /// sitting a hair above its bound and stepping outward throttles the step
    /// for ALL variables (a documented hard-stall mechanism in GEMS3K's Optima
    /// integration - a trace species 2e-19 above a 1e-13 floor produced
    /// betamin ~ 8e-12 and froze the whole iterate).
    ///
    /// When set to a value `f > 0`, ORCHESTRA's two-part rule applies instead:
    ///   1. a variable whose OWN beta is <= `f` is EXCLUDED from setting the
    ///      shared `betamin` (it no longer drags everyone else down), and
    ///   2. that variable is STILL clamped to its own (smaller) beta when the
    ///      step is applied - it moves conservatively on its own scale rather
    ///      than at the now-larger shared factor, so it can never blow through
    ///      its own bound.
    /// Equivalently, every variable moves by `min(own_beta, betamin)`.
    ///
    /// Note this makes the applied step a per-variable rescaling rather than a
    /// pure shrink along the Newton ray - the direction is no longer exactly
    /// preserved. ORCHESTRA accepts the same trade deliberately (see its own
    /// in-code comment "is it better to update unknowns that are very
    /// sensitive with their own small factor?"). The duals `w` have no bounds
    /// and always use the shared `betamin`.
    ///
    /// MEASURED AND REJECTED as a default, 2026-08-26 - kept only as
    /// default-off infrastructure so the experiment is cheap to re-run.
    /// The pathology it targets is real and pervasive (on GEMS3K's
    /// f_GEOTHERM, 764 of 1765 AOP iterations had a shared betamin below
    /// 1e-6, i.e. the step was annihilated by a single trace outlier), and
    /// the rule does help the stationarity residual - `||ex||max` decays
    /// slightly FASTER with it on. But it breaks the linear equality
    /// constraint `Aex*x = be`: a uniform blend keeps `Aex*x` a convex
    /// combination of `Aex*xo` and `be`, so feasibility is preserved,
    /// whereas a per-variable rescaling lands on no such line. Measured
    /// directly on GEMS3K's f_Solvus at f = 1e-9: `||ew||max` jumps from
    /// ~105 to ~10000 by iteration 5 and never recovers, and `Error`
    /// becomes dominated by the constraint residual rather than by
    /// optimality. Across a 17-project sweep at f = 1e-5, 9 projects
    /// improved by up to 34% while 5 went from converged to failed with
    /// GEMS3K's own mass-balance check firing; the response is
    /// non-monotonic in `f` (one project ran OK/145 -> OK/102 -> FAIL ->
    /// OK/204 -> OK/182 over f = 0, 1e-12, 1e-9, 1e-7, 1e-5), i.e. the
    /// wins are "the next Newton step happened to repair feasibility in
    /// time", not a sound mechanism.
    ///
    /// This is an architectural mismatch, not a tuning problem: ORCHESTRA's
    /// unknowns each own a separate residual equation and are not coupled by
    /// a shared linear equality, so it can rescale them independently at no
    /// cost. Any per-variable step-length scheme in this formulation has the
    /// same defect - which also retroactively explains the earlier, simpler
    /// variant tried on the GEMS3K side ("exclude a tiny beta from betamin
    /// and let the final box clamp handle that variable"), whose regression
    /// was recorded but not explained at the time: the clamp is per-variable
    /// too, so it breaks feasibility by exactly the same mechanism.
    double min_common_beta = 0.0;
};

} // namespace Optima
