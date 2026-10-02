# Optima options and methods — a short guide

Optima solves constrained optimization problems: minimize `f(x, p)` subject to linear and
nonlinear equality and inequality constraints, with lower and upper bounds on `x` and `p`.
This guide lists what you can set and when each setting matters.

## Using the solver

| Method | What it does |
|---|---|
| `Solver::setOptions(options)` | Sets the options used by the next solves. |
| `Solver::solve(problem, state)` | Solves the problem, starting from `state`, and writes the solution back into `state`. |
| `Solver::solve(problem, state, sensitivity)` | The same, plus derivatives of the solution with respect to the parameters `c`. |

**Starting point.** `state` is both input and output. A good start (for example, the solution of
a nearby problem) usually means fewer iterations. Solving a series of similar problems by reusing
the previous state is the main way to gain speed.

## The problem — `Problem`

| Field | What it is |
|---|---|
| `f` | Objective function: value, gradient and Hessian. |
| `he`, `hg` | Nonlinear equality (`he = 0`) and inequality (`hg >= 0`) constraints. |
| `v` | External nonlinear constraints `v(x, p) = 0` (used for the parameters `p`). |
| `r` | Optional function that computes data shared by the objective and constraints once per point. |
| `be`, `bg` | Right-hand sides of the linear equality and inequality constraints. |
| `xlower`, `xupper`, `plower`, `pupper` | Bounds on `x` and `p`. |
| `c`, `bec`, `bgc` | Parameters and derivatives of `be`, `bg` with respect to them, for sensitivities. |

## The answer — `State` and `Result`

`State` holds the solution `x`, `p`, the Lagrange multipliers (`ye`, `yg`, `ze`, `zg`) and the
stability measures `s`. A variable at its lower bound whose `s` is positive would not want to
grow, so it is correctly at the bound. `State` also lists stable (`js`) and unstable (`ju`) variables.

`Result` says whether the solve succeeded (`succeeded`, `failure_reason`) and reports the number of
iterations, the final errors, and counts and timings of function evaluations. Check `succeeded`
before using the answer.

## Options — `Options`

### General

| Option | Default | When it matters |
|---|---|---|
| `maxiters` | 200 | Hard problems or poor starting points may need more iterations. Raise it if solves stop with "Max iterations reached". |
| `convergence.tolerance` | 1e-8 | The solve stops when the error falls below this. Smaller = more precise and more iterations. Match it to the scale of your problem. |
| `convergence.check` | none | An extra test of your own. The solve also stops, as converged, when it returns true. Useful when your problem has its own natural stopping test. |
| `newtonstep.linearsolver.method` | `Nullspace` | How each step's linear system is solved. `Fullspace`: small problems. `Nullspace`: general default, good for a dense Hessian with many constraints. `Rangespace`: only when the Hessian is diagonal. |
| `output.active`, `output.filename`, ... | off | Writes a per-iteration table, for debugging. `xnames`, `pnames`, ... label the columns. |

### Keeping steps inside the bounds — `backtracksearch`

Each step is shortened so that no variable crosses its bound. One common factor is used for all
variables, so the step direction is kept.

| Option | Default | When it matters |
|---|---|---|
| `apply_min_max_fix_and_accept` | false | Instead of shortening the step, clip each variable to its bounds. Can help when one variable near its bound keeps the step tiny. The step direction is changed. |

### Line search — `linesearch`

When a step makes the error more than `trigger...by_factor` times larger than before (default:
more than double), the line search looks at all points between the old point and the new one and
keeps the one with the smallest error. This is usually a shorter step, but it can be the full
step. The "new point" is the step after it was already cut to stay within the bounds.
Off by default.

| Option | Default | When it matters |
|---|---|---|
| `enabled` | false | Turn on for problems where full steps overshoot and the solve oscillates or diverges. |
| `use_unmasked_error` | false | Judge progress by the full error, including variables sitting on a bound. Recommended when the line search is on and many variables are at bounds. |
| `trigger_when_current_error_is_greater_than_previous_error_by_factor` | 2.0 | How much worse the error must get before a line search runs. Lower = more line searches (safer, more cost). |
| `tolerance`, `maxiterations` | 1e-5, 20 | Precision and effort of each line search. |
| `stall_escape_after` | 0 (off) | If several line searches in a row get nowhere, take one full step to break out. Helps when the solve freezes at a constant error. |
| `stall_escape_tolerance` | 1e-8 | How small a change counts as "getting nowhere". |
| `nonmonotone_window` | 0 (off) | Allow the error to rise for a few steps before it counts as worse. Can speed up slow, zig-zagging solves; can also make others less reliable. |
| `reject_if_worse` | false | If the shorter step did not help either, use the full step. Helps when the line search keeps making tiny steps in the wrong direction. |

Line-search options can be tested from the code with the macros `OPTIMA_LINESEARCH_STALL_ESCAPE`
and `OPTIMA_LINESEARCH_REJECT_WORSE`.

### Options that are currently not read by the solver

`errorstatus.*`, `steepestdescent.*` and
`linesearch.trigger_when_current_error_is_greater_than_initial_error_by_factor` exist but do not
change the result in this version.

## Diagnostics

- Set the environment variable `OPTIMA_BETA_PROBE` to a file name to log, for every step, the
  common step factor and the variable that limited it. Use it when a solve makes almost no
  progress per iteration.
- When the library loads, it prints `[Optima] LOCAL MODIFIED BUILD loaded` to stderr, so you can
  confirm this build, and not a stock Optima, is in use.

## Quick choices

| Situation | Try |
|---|---|
| Stops at the iteration limit, error still falling | Raise `maxiters`; give a better starting state. |
| Error oscillates or grows | `linesearch.enabled = true` with `use_unmasked_error = true`. |
| Line search on, error frozen for many iterations | `stall_escape_after = 10`. |
| Line search on, many tiny steps the wrong way | `reject_if_worse = true`. |
| Almost no progress per step, many variables at bounds | Check with `OPTIMA_BETA_PROBE`; try `apply_min_max_fix_and_accept`. |
