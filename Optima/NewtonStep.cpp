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

#include "NewtonStep.hpp"

// Optima includes
#include <Optima/Exception.hpp>
#include <Optima/LinearSolver.hpp>
#include <Optima/Trace.hpp>

#include <cmath>

namespace Optima {

struct NewtonStep::Impl
{
    MasterDims dims;           ///< The dimensions of the master variables.
    LinearSolver linearsolver; ///< The linear solver for the master matrix equations.
    MasterVector du;           ///< The Newton step for master variables u = (x, p, w).
    LinearSolverMethod method = LinearSolverMethod::Nullspace; ///< The linear solver method in use.
    Vector xlower;             ///< The lower bounds for variables *x*.
    Vector xupper;             ///< The upper bounds for variables *x*.
    Vector plower;             ///< The lower bounds for variables *p*.
    Vector pupper;             ///< The upper bounds for variables *p*.

    Impl()
    {}

    auto setOptions(const NewtonStepOptions options) -> void
    {
        linearsolver.setOptions(options.linearsolver);
        method = options.linearsolver.method;
    }

    auto initialize(const MasterProblem& problem) -> void
    {
        dims = problem.dims;
        xlower = problem.xlower;
        xupper = problem.xupper;
        plower = problem.plower;
        pupper = problem.pupper;
        du.resize(dims);
    }

    auto apply(const ResidualFunction& F, MasterVectorView uo, MasterVectorRef u) -> void
    {
        sanitycheck();
        const auto res = F.result();
        const auto Jc = res.Jc;
        const auto Fc = res.Fc;
        linearsolver.decompose(Jc);
        linearsolver.solve(Jc, Fc, du);

        if(auto* fp = trace::file())
        {
            // The Newton step: method and the largest change in x, p and w.
            Index idu = 0;
            double nx = 0.0;
            for(Index i = 0; i < du.x.size(); ++i)
                if(std::abs(du.x[i]) > nx) { nx = std::abs(du.x[i]); idu = i; }
            const double np = du.p.size() ? du.p.cwiseAbs().maxCoeff() : 0.0;
            const double nw = du.w.size() ? du.w.cwiseAbs().maxCoeff() : 0.0;
            const char* m = method == LinearSolverMethod::Nullspace ? "nullspace" : method == LinearSolverMethod::Rangespace ? "rangespace" : "fullspace";
            fprintf(fp, "NEWTON solve=%ld attempt=%d it=%ld method=%s ndu_x=%.6e ndu_p=%.6e ndu_w=%.6e idu=%ld\n",
                trace::state().solve, trace::state().attempt, trace::state().iter, m, nx, np, nw, (long)idu);
            fflush(fp);
        }
        u.x.noalias() = uo.x + du.x;
        u.p.noalias() = uo.p + du.p;
        u.w.noalias() = uo.w + du.w;
    }

    auto sanitycheck() const -> void
    {
        assert(xlower.size() == dims.nx);
        assert(xupper.size() == dims.nx);
        assert(plower.size() == dims.np);
        assert(pupper.size() == dims.np);
    }
};

NewtonStep::NewtonStep()
: pimpl(new Impl())
{}

NewtonStep::NewtonStep(const NewtonStep& other)
: pimpl(new Impl(*other.pimpl))
{}

NewtonStep::~NewtonStep()
{}

auto NewtonStep::operator=(NewtonStep other) -> NewtonStep&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto NewtonStep::setOptions(const NewtonStepOptions& options) -> void
{
    pimpl->setOptions(options);
}

auto NewtonStep::initialize(const MasterProblem& problem) -> void
{
    pimpl->initialize(problem);
}

auto NewtonStep::apply(const ResidualFunction& F, MasterVectorView uo, MasterVectorRef u) -> void
{
    pimpl->apply(F, uo, u);
}

} // namespace Optima
