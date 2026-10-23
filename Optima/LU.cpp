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

#include "LU.hpp"

// C++ includes
#include <cassert>
#include <cmath>
#include <limits>

// Eigen includes
#include <Eigen/src/LU/FullPivLU.h>
#include <Eigen/src/LU/PartialPivLU.h>

// Optima includes
#include <Optima/Macros.hpp>
#include <Optima/Trace.hpp>

namespace Optima {

struct LU::Impl
{
    //======================================================================
    // Note: The full pivoting strategy is needed at the moment to resolve
    // singular matrices. Using a partial pivoting scheme via PartialPivLU
    // would need to be combined with a search for linearly dependent rows in
    // the produced upper triangular matrix U.
    //======================================================================

    /// The base LU solver from Eigen library.
    Eigen::FullPivLU<Matrix> lu;

    /// The partial-pivot LU, used for a matrix of at least partialMinSize rows.
    Eigen::PartialPivLU<Matrix> pplu;

    /// The matrix last given to decompose, kept while the partial-pivot LU is in use (for the fallback).
    Matrix Acopy;

    /// The smallest matrix factorised with partial pivoting (0 = never).
    Index partialMinSize = 0;

    /// True while the current factorisation is the partial-pivot one.
    bool usePartial = false;

    /// The workspace for y in the partial-pivot solve.
    Vector yp;

    /// The workspace for modified matrix U.
    Matrix U;

    /// The flags that indicate if an equation is linearly independent (non-zero value).
    Indices is_li;

    /// The rank of the lineary system, not of the coefficient matrix (depends on right-hand side vector!)
    Index rank;

    /// Construct a default Impl object.
    Impl()
    {}

    /// Construct an Impl object with given matrix.
    Impl(MatrixView A)
    {
        decompose(A);
    }

    /// Return true if empty.
    auto empty() const -> bool
    {
        return usePartial ? pplu.matrixLU().size() == 0 : lu.matrixLU().size() == 0;
    }

    /// Compute the LU decomposition of the given matrix.
    auto decompose(MatrixView A) -> void
    {
        const auto m = A.rows();
        const auto n = A.cols();
        assert(n == m);
        usePartial = partialMinSize > 0 && n >= partialMinSize;
        if(usePartial)
        {
            Acopy = A;
            pplu.compute(Acopy);
        }
        else lu.compute(A);
    }

    /// Solve the linear system `Ax = b` using the LU decomposition obtained with @ref decompose.
    auto solve(VectorView b, VectorRef x) -> void
    {
        x = b;
        solve(x);
    }

    /// Solve the linear system `Ax = b` using the LU decomposition obtained with @ref decompose.
    auto solve(VectorRef x) -> void
    {
        if(usePartial && solvePartial(x))
            return;
        const auto& n  = lu.rows();
        const auto& P  = lu.permutationP();
        const auto& Q  = lu.permutationQ();
        const auto& M  = lu.matrixLU();
        const auto& Lv = M.triangularView<Eigen::UnitLower>();
        const auto& Uv = U.triangularView<Eigen::Upper>();

        assert(n == x.rows());

        P.applyThisOnTheLeft(x);
        x = Lv.solve(x);
        assembleU(x);
        x = Uv.solve(x);
        Q.applyThisOnTheLeft(x);
        Q.applyThisOnTheLeft(is_li);

        // TODO; In LU, x should have +inf or -inf to indicate extremely large steps and their directions. Then a line search would be used to find a reasonable step length/
    }

    /// Solve with the partial-pivot factors. Returns false, with x unchanged and the full-pivot factors computed,
    /// when assembleU's test would drop an equation (|U_ii| <= eps |y_i|): that step goes the full-pivot way.
    auto solvePartial(VectorRef x) -> bool
    {
        const auto n = pplu.rows();
        const auto& M = pplu.matrixLU();
        const auto eps = std::numeric_limits<double>::epsilon();
        yp.noalias() = pplu.permutationP() * x;
        M.triangularView<Eigen::UnitLower>().solveInPlace(yp);
        trace::state().lupartial++;
        for(Index i = 0; i < n; ++i)
        {
            if(std::abs(M(i, i)) <= eps * std::abs(yp[i]))
            {
                trace::state().lufallbacks++;
                lu.compute(Acopy);
                usePartial = false; // later solves with this matrix use the full-pivot factors too
                return false;
            }
        }
        M.triangularView<Eigen::Upper>().solveInPlace(yp);
        x = yp;
        is_li.setOnes(n);
        rank = n;
        return true;
    }

    /// Assemble the U matrix with given y, where y is the solution of L*y = P*b.
    auto assembleU(VectorRef y) -> void
    {
        const auto n = lu.rows();

        U.resize(n, n);
        U = lu.matrixLU();

        is_li.setOnes(n); // set all equations as linearly independent to start with

        // Skip the rest if there is only one equation.
        if(n == 1)
            return;

        const auto D = U.diagonal().cwiseAbs();
        const auto eps = std::numeric_limits<double>::epsilon();

        using std::max;
        using std::abs;

        rank = n; // start full rank, decrease as we go along through the diagonal of U (from the bottom!)
        for(auto i = 1; i <= n; ++i)
        {
            // Check diagonal entry is not very small compared to other entries
            // on same row of U and also right-hand side entry in y. The idea
            // is that if a diagonal entry is very small, but the other
            // coefficients along the same row of U is also very small and the
            // corresponding entry in y is equally small, then it is safe not
            // to discard this linear equation. Otherwise, we discard it, to
            // avoid extremely large values when we divide a larger number by
            // the diagonal pivot (very small).
            if(D[n - i] <= eps * abs(y[n - i]))
            {
                U.row(n - i).tail(n - i).fill(0.0); // avoid going to the lower triangular part
                U.col(n - i).head(n - i).fill(0.0);
                U(n - i, n - i) = 1.0; // diagonal entry is 1 while the rest entries along the same row/col are 0
                y[n - i] = 0.0; // the solution of the corresponding unknown should be zero
                is_li[n - i] = 0; // equation not linearly independent
                --rank;
            }
        }
    }
};

LU::LU()
: pimpl(new Impl())
{}

LU::LU(const LU& other)
: pimpl(new Impl(*other.pimpl))
{}

LU::~LU()
{}

auto LU::operator=(LU other) -> LU&
{
    pimpl = std::move(other.pimpl);
    return *this;
}

auto LU::empty() const -> bool
{
    return pimpl->empty();
}

auto LU::setPartialMinSize(Index n) -> void
{
    pimpl->partialMinSize = n;
}

auto LU::decompose(MatrixView A) -> void
{
    pimpl->decompose(A);
}

auto LU::solve(VectorView b, VectorRef x) -> void
{
    x = b;
    pimpl->solve(x);
}

auto LU::solve(VectorRef x) -> void
{
    pimpl->solve(x);
}

auto LU::rank() const -> Index
{
    return pimpl->rank;
}

auto LU::matrixLU() const -> MatrixView
{
    return pimpl->usePartial ? MatrixView(pimpl->pplu.matrixLU()) : MatrixView(pimpl->lu.matrixLU());
}

auto LU::P() const -> PermutationMatrix
{
    return pimpl->usePartial ? PermutationMatrix(pimpl->pplu.permutationP()) : PermutationMatrix(pimpl->lu.permutationP());
}

auto LU::Q() const -> PermutationMatrix
{
    // The partial-pivot factors permute rows only.
    if(pimpl->usePartial) { PermutationMatrix q(pimpl->pplu.rows()); q.setIdentity(); return q; }
    return pimpl->lu.permutationQ();
}

} // namespace Optima
