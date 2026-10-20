# Optima options and methods — a short guide

Optima solves constrained optimization problems. It minimizes `f(x, p)` with linear and nonlinear equality
and inequality constraints, and with lower and upper bounds on `x` and `p`. This guide tells you which
settings you can change, and when each setting is important.

## Using the solver

| Method | What it does |
|---|---|
| `Solver::setOptions(options)` | Sets the options for the next solves. |
| `Solver::solve(problem, state)` | Solves the problem from `state` and writes the solution into `state`. |
| `Solver::solve(problem, state, sensitivity)` | Does the same and also calculates the derivatives of the solution with respect to the parameters `c`. |

**Starting point.** `state` is the input and the output. A good start, for example the solution of a
similar problem, usually needs fewer iterations. To solve a series of similar problems quickly, use the
previous state as the start of the next solve.

## The problem — `Problem`

| Field | What it is |
|---|---|
| `f` | The objective function: value, gradient and Hessian. |
| `he`, `hg` | Nonlinear equality constraints (`he = 0`) and inequality constraints (`hg >= 0`). |
| `v` | External nonlinear constraints `v(x, p) = 0`, for the parameters `p`. |
| `r` | An optional function. It calculates, once for each point, the data that the objective and the constraints share. |
| `be`, `bg` | The right-hand sides of the linear equality and inequality constraints. |
| `xlower`, `xupper`, `plower`, `pupper` | The bounds on `x` and `p`. |
| `c`, `bec`, `bgc` | The parameters, and the derivatives of `be` and `bg` with respect to them, for sensitivities. |

## The answer — `State` and `Result`

`State` holds the solution `x` and `p`, the Lagrange multipliers (`ye`, `yg`, `ze`, `zg`) and the stability
measures `s`. If a variable is at its lower bound and its `s` is positive, the variable does not tend to grow.
Thus it is correctly at the bound. `State` also lists the stable variables (`js`) and the unstable variables (`ju`).

`Result` tells you whether the solve succeeded (`succeeded`, `failure_reason`). It also gives the number of
iterations, the final errors, and the counts and times of the function evaluations. Check `succeeded` before
you use the answer.

## Options — `Options`

### General

| Option | Default | When it is important |
|---|---|---|
| `maxiters` | 200 | A difficult problem or a poor start can need more iterations. If solves stop with "Max iterations reached", increase it. |
| `convergence.tolerance` | 1e-8 | The solve stops when the error is less than this value. A smaller value gives a more precise answer and more iterations. Use a value that agrees with the scale of your problem. |
| `convergence.check` | none | Your own stop test. When it returns true, the solve stops and counts as converged. Use it if your problem has its own natural stop test. |
| `newtonstep.linearsolver.method` | `Nullspace` | The method that solves the linear system of each step. `Fullspace`: for small problems. `Nullspace`: the general method, good for a dense Hessian with many constraints. `Rangespace`: only for a diagonal Hessian. |
| `newtonstep.linearsolver.lu_partial_min_size` | 0 (off) | `Nullspace` only. The solver factorizes a matrix of this size or larger with partial pivoting, which is faster. If a step would lose an equation, the solver does that step again with full pivoting. 0 means full pivoting always. Use about 100 for large systems. |
| `output.active`, `output.filename`, ... | off | Writes a table with one line for each iteration, for debugging. `xnames`, `pnames`, ... give the column labels. |

### Keeping steps inside the bounds — `backtracksearch`

Each step is made shorter so that no variable goes across its bound. All variables use the same factor, so
the direction of the step does not change.

| Option | Default | When it is important |
|---|---|---|
| `apply_min_max_fix_and_accept` | false | Moves each variable that is outside its bounds back to the bound, instead of a shorter step. This can help when one variable near its bound keeps the step very small. This option changes the direction of the step. |
| `snap_leaving_after` | 0 (off) | A variable can limit the step many times in a row on its way to its lower bound. Its target is then far below the bound. After this number of steps, the solver puts the variable on the bound. The other variables then limit the step. Use 2 for warm starts, where a disappearing species can otherwise go down slowly. |

### Line search — `linesearch`

The line search is off by default. When it is on, it starts if a step makes the error too large. The default
limit is 2 times the previous error. It then looks at all the points between the old point
and the new point and keeps the point with the smallest error. Usually this is a shorter step, but it can be
the full step. The new point is the step after it was made shorter to stay inside the bounds.

| Option | Default | When it is important |
|---|---|---|
| `enabled` | false | Set it to true for problems where full steps go too far and the solve oscillates or diverges. |
| `use_unmasked_error` | false | Uses the full error, which includes the variables on a bound. Use it when the line search is on and many variables are on bounds. |
| `trigger_when_current_error_is_greater_than_previous_error_by_factor` | 2.0 | How much larger the error must become before a line search starts. A smaller value gives more line searches: safer, but slower. |
| `tolerance`, `maxiterations` | 1e-5, 20 | The precision and the maximum work of each line search. |
| `stall_escape_after` | 0 (off) | After this number of line searches in a row that do not change the error, the solve uses one full step. Use it when the error stays the same for many iterations. |
| `stall_escape_tolerance` | 1e-8 | The relative change of the error that counts as no change. |
| `reject_if_worse` | false | If a line search does not make the error less than before the step, the solve uses the full step. Use it when the line search makes many very small steps in the wrong direction. |

Client code can test for the line-search options with the macros `OPTIMA_LINESEARCH_STALL_ESCAPE` and
`OPTIMA_LINESEARCH_REJECT_WORSE`.

### Options that the solver does not read now

`errorstatus.*`, `steepestdescent.*` and
`linesearch.trigger_when_current_error_is_greater_than_initial_error_by_factor` exist, but they do not change
the result in this version.

## Selecting the basic variables

At each step the solver selects one basic variable for each independent constraint. Two tests keep this
selection independent. They are on by default and have no option. Use the environment variables only for tests.

| Environment variable | Default | What it does |
|---|---|---|
| `OPTIMA_ECH_RELPIVOT` | 1e-8 | A candidate for a swap must be larger than this fraction of the largest entry in its row. 0 turns this test off. |
| `OPTIMA_ECH_RANKCHECK` | on | A swap is made only if the new basic columns keep full rank. 0 turns this test off. If a swap fails the test, the solver tries the next candidate. |
| `OPTIMA_ECH_RANKTOL` | 1e-8 | For the rank test, a pivot smaller than this fraction of the largest pivot counts as zero. |

## Diagnostics

Each switch is an environment variable that gives a file name. It writes only when you set it.

| Environment variable | What it writes |
|---|---|
| `OPTIMA_BETA_PROBE` | For each step, the final common step factor and the variable that limited it. Use it when a solve makes almost no progress in each iteration. |
| `OPTIMA_LU_PROBE` | For each linear solve, a comparison of full and partial pivoting: pivot range, dropped equations, residuals and time. |
| `GEMS3K_OPTIMA_STEPTRACE_FILE` | A debug record of each iteration: the Newton step size, the step factor and the variable that limited it, the objective and the errors, and the reason the solve stopped. |

When the library loads, it writes `[Optima] LOCAL MODIFIED BUILD loaded` to stderr. Thus you can make sure
that your program uses this build and not a standard Optima.

## Quick choices

| Problem | Try this |
|---|---|
| The solve stops at the iteration limit and the error still decreases. | Increase `maxiters`. Give a better start. |
| The error oscillates or increases. | Set `linesearch.enabled = true` and `use_unmasked_error = true`. |
| The line search is on and the error does not change for many iterations. | Set `stall_escape_after = 10`. |
| The line search is on and makes many small steps in the wrong direction. | Set `reject_if_worse = true`. |
| Each step makes almost no progress and many variables are on bounds. | Look at `OPTIMA_BETA_PROBE`. Try `apply_min_max_fix_and_accept`. |
| A warm start is slow because a species goes down to its bound a little at each step. | Set `snap_leaving_after = 2`. |
| Large problems are slow in the linear solve. | Set `lu_partial_min_size = 100`. |
