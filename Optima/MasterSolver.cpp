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

#include "MasterSolver.hpp"

// Optima includes
#include <Optima/Convergence.hpp>
#include <Optima/ErrorControl.hpp>
#include <Optima/Exception.hpp>
#include <Optima/MasterProblem.hpp>
#include <Optima/MasterVector.hpp>
#include <Optima/NewtonStep.hpp>
#include <Optima/Options.hpp>
#include <Optima/Outputter.hpp>
#include <Optima/ResidualErrors.hpp>
#include <Optima/ResidualFunction.hpp>
#include <Optima/Result.hpp>
#include <Optima/SensitivitySolver.hpp>
#include <Optima/TransformStep.hpp>
#include <Optima/Trace.hpp>

#include <algorithm>

namespace Optima {

const auto CONTINUE = true;
const auto STOP     = false;

struct MasterSolver::Impl
{
    MasterDims dims;
    ResidualFunction F;
    ResidualErrors E;
    MasterVector uo;
    NewtonStep newtonstep;
    TransformStep transformstep;
    ErrorControl errorcontrol;
    Convergence convergence;
    bool converged = false;
    SensitivitySolver sensitivitysolver;
    Outputter outputter; ///< The object used to output the current state of the computation.
    Result result;
    Options options;
    trace::Recorder tracer;             ///< The ITER and STOP records, written only while tracing.

    Impl()
    {}

    /// Write one ITER record for the current state.
    auto traceIter() -> void
    {
        if(!trace::file()) return;
        const auto& r = F.result();
        const auto& jb = r.Jm.RWQ.jb;
        const auto& jn = r.Jm.RWQ.jn;
        const auto& ss = r.stabilitystatus;
        tracer.iter((long)result.iterations, (double)r.f.f, E.error(), E.errorx(), E.errorp(), E.errorw(),
            jb.data(), (long)jb.size(), jn.data(), (long)jn.size(), uo.x.data(),
            (long)r.Jc.dims.nbs, (long)r.Jc.dims.nl, (long)ss.js.size(), (long)ss.jlu.size(), (long)ss.juu.size());
    }

    /// Write the STOP record: why the loop ended.
    auto traceStop() -> void
    {
        if(!trace::file()) return;
        const char* reason = converged ? (E.error() < options.convergence.tolerance ? "tol" : "check")
                                       : (result.iterations > options.maxiters ? "maxiters" : "fatal");
        tracer.stop((long)result.iterations, reason, E.error(), options.convergence.tolerance, options.maxiters,
            static_cast<bool>(options.convergence.check));
    }

    auto outputHeaderTop() -> void
    {
        if(!options.output.active) return;

        outputter.addEntry("Iteration");
        outputter.addEntry("f");
        outputter.addEntry("Error");
        outputter.addEntry("||ex||max");
        outputter.addEntry("||ep||max");
        outputter.addEntry("||ew||max");
        outputter.addEntries("x", dims.nx, options.output.xnames);
        outputter.addEntries("p", dims.np, options.output.pnames);
        outputter.addEntries("y", dims.ny, options.output.ynames);
        outputter.addEntries("z", dims.nz, options.output.znames);
        outputter.addEntries("s", dims.nx, options.output.xnames);
        outputter.addEntries("ex", dims.nx, options.output.xnames);
        outputter.addEntries("ep", dims.np, options.output.pnames);
        outputter.addEntries("ey", dims.ny, options.output.ynames);
        outputter.addEntries("ez", dims.nz, options.output.znames);
        outputter.addEntry("Basic Variables");
        outputter.outputHeader();
    };

    auto outputHeaderBottom() -> void
    {
        if(!options.output.active) return;
        outputter.outputHeader();
    };

    auto outputCurrentState() -> void
    {
        if(!options.output.active) return;
        const auto& Fres = F.result();
        const auto& jb = Fres.Jm.RWQ.jb;
        const auto& xnames = options.output.xnames;
        outputter.addValue(result.iterations);
        outputter.addValue(Fres.f.f);
        outputter.addValue(E.error());
        outputter.addValue(E.errorx());
        outputter.addValue(E.errorp());
        outputter.addValue(E.errorw());
        outputter.addValues(uo.x);
        outputter.addValues(uo.p);
        outputter.addValues(uo.w.head(dims.ny));
        outputter.addValues(uo.w.tail(dims.nz));
        outputter.addValues(Fres.stabilitystatus.s);
        outputter.addValues(E.ex());
        outputter.addValues(E.ep());
        outputter.addValues(E.ew().head(dims.ny));
        outputter.addValues(E.ew().tail(dims.nz));
        std::stringstream ss;
        if(xnames.empty()) for(auto i : jb) ss << i << " ";
        else for(auto i : jb) ss << xnames[i] << " ";
        outputter.addValue(ss.str());
        outputter.outputState();
    };

    auto solve(const MasterProblem& problem, MasterState& state) -> Result
    {
        auto& u = state.u;
        initialize(problem, u);
        do step(u); while(stepping(u));
        finalize(state);
        return result;
    }

    auto solve(const MasterProblem& problem, MasterState& state, MasterSensitivity& sensitity) -> Result
    {
        solve(problem, state);
        F.updateOnlyJacobian(state.u); // update the Jacobian matrices wrt x, p, c
        sensitivitysolver.solve(F, state, sensitity);
        return result;
    }

    auto setOptions(const Options& opts) -> void
    {
        options = opts;
        newtonstep.setOptions(opts.newtonstep);
        convergence.setOptions(opts.convergence);
        errorcontrol.setOptions({opts.errorstatus, opts.backtracksearch, opts.linesearch});
        outputter.setOptions(opts.output);
    }

    auto initialize(const MasterProblem& problem, MasterVectorRef u) -> void
    {
        sanitycheck(problem, u);
        dims = problem.dims;
        result = {};
        if(trace::file())
        {
            trace::state().iter = 0;
            if(trace::state().auto_attempt) trace::state().attempt = ++trace::state().solves_in_call;
            tracer.begin(problem.xlower.data(), problem.xupper.data(), (long)problem.xlower.size());
        }
        u.x.noalias() = min(max(u.x, problem.xlower), problem.xupper);
        exXlo = problem.xlower; exXup = problem.xupper; exD0.resize(0); exR0 = -1.; exC0 = 0.;
        u.p.noalias() = min(max(u.p, problem.plower), problem.pupper);
        uo = u;
        F.initialize(problem);
        F.update(u);
        E.initialize(problem);
        E.update(u, F);
        transformstep.initialize(problem);
        newtonstep.initialize(problem);
        errorcontrol.initialize(problem);
        convergence.initialize(problem);
        sensitivitysolver.initialize(problem);
        outputter.clear();
        outputHeaderTop();
    }

    auto stepping(MasterVectorRef u) -> bool
    {
        convergence.update(E);
        // Convergence is checked before the iteration limit, so the final step that
        // the limit allows (iterations == maxiters + 1) is reported as converged if it is.
        ConvergenceCheckArgs args{dims, F, E, uo, u, result};
        converged = convergence.converged(args);
        uo = u;
        if(converged)
            return STOP;
        if(result.iterations > options.maxiters)
            return STOP;
        return CONTINUE;
    }

    auto step(MasterVectorRef u) -> void
    {
        outputCurrentState();
        traceIter();
        result.iterations += 1;
        trace::state().iter = result.iterations;
        newtonstep.apply(F, uo, u);
        transformstep.execute(uo, u, F, E);
        errorcontrol.execute(uo, u, F, E);
        F.update(u);
        E.update(u, F);
        extrapolate(u);
    }

    // Jump ahead (options.jumpahead): when two successive steps point the same way (cosine > 0.99) and shrink
    // at a steady ratio r in [0.5, 1), try x + (g - 1) d with g = 1/(1 - r), g <= 1000. The jump is skipped if
    // any variable would cross a bound, and it is kept only if the error falls.
    Vector exD0, exXlo, exXup;
    double exR0 = -1., exC0 = 0.;
    auto extrapolate(MasterVectorRef u) -> void
    {
        if(!options.jumpahead) return;
        const Vector d = u.x - uo.x;
        const double nd = d.norm();
        double r = -1., cs = 0.;
        if(exD0.size() == d.size() && exD0.norm() > 0. && nd > 0.)
        {
            r = nd / exD0.norm();
            cs = d.dot(exD0) / (nd * exD0.norm());
        }
        bool fired = false;
        if(cs > 0.99 && exC0 > 0.99 && r >= 0.5 && r < 1.0 && exR0 >= 0.5 && exR0 < 1.0001)
        {
            const double g = 1. / (1. - std::min({r, exR0, 0.9999}));
            bool inside = g <= 1000.;
            for(Index i = 0; i < d.size() && inside; i++)
            {
                const double xn = u.x[i] + (g - 1.) * d[i];
                if(xn < exXlo[i] || xn > exXup[i]) inside = false;
            }
            if(inside)
            {
                const double e0 = E.error();
                const Vector xkeep = u.x;
                u.x.noalias() = u.x + (g - 1.) * d;
                F.update(u);
                E.update(u, F);
                if(E.error() < e0) fired = true;
                else { u.x = xkeep; F.update(u); E.update(u, F); }
            }
        }
        // After a jump the history starts again.
        if(fired) { exD0.resize(0); exR0 = -1.; exC0 = 0.; }
        else { exD0 = d; exR0 = r; exC0 = cs; }
    }

    auto finalize(MasterState& state) -> void
    {
        // Set the final values of the result object
        result.succeeded = converged;
        result.failure_reason = "";
        result.failure_fatal = false;
        result.error = E.error();
        result.error_optimality = E.errorx();
        result.error_feasibility = E.errorw();

        if(!converged) {
            if(result.iterations > options.maxiters) {
                result.failure_reason = "Max iterations reached.";
            }
            else {
                result.failure_reason = "An unknown fatal error occurred.";
                result.failure_fatal = true;
            }
        }

        outputCurrentState();
        traceIter();
        traceStop();
        outputHeaderBottom();
        auto const& Fresult = F.result();
        auto const& ss = Fresult.stabilitystatus;
        state.s = ss.s;
        state.js = ss.js;
        state.ju = ss.ju;
        state.jlu = ss.jlu;
        state.juu = ss.juu;
        state.jb = Fresult.Jc.jb;
        state.jn = Fresult.Jc.jn;
    }

    auto sanitycheck(const MasterProblem& problem, MasterVectorRef u) -> void
    {
        assert(problem.f.initialized());
        assert(problem.dims.nz == 0 || problem.h.initialized());
        assert(problem.dims.np == 0 || problem.v.initialized());
        assert(problem.dims.nx == problem.xlower.size());
        assert(problem.dims.nx == problem.xupper.size());
        assert(problem.dims.np == problem.plower.size());
        assert(problem.dims.np == problem.pupper.size());
        assert(problem.dims.nx == u.x.size());
        assert(problem.dims.np == u.p.size());
        assert(problem.dims.nw == u.w.size());
    }
};

MasterSolver::MasterSolver()
: pimpl(new Impl())
{}

MasterSolver::MasterSolver(const MasterSolver& other)
: pimpl(new Impl(*other.pimpl))
{}

MasterSolver::~MasterSolver()
{}

auto MasterSolver::operator=(MasterSolver other) -> MasterSolver&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto MasterSolver::setOptions(const Options& options) -> void
{
    pimpl->setOptions(options);
}

auto MasterSolver::solve(const MasterProblem& problem, MasterState& state) -> Result
{
    return pimpl->solve(problem, state);
}

auto MasterSolver::solve(const MasterProblem& problem, MasterState& state, MasterSensitivity& sensitivity) -> Result
{
    return pimpl->solve(problem, state, sensitivity);
}

} // namespace Optima
