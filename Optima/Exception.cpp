// Optima is a C++ library for solving linear and non-linear constrained optimization problems.
//
// Copyright © 2020-2024 Allan Leal
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library. If not, see <http://www.gnu.org/licenses/>.

#include "Exception.hpp"

// C++ includes
#include <memory>
#include <mutex>

namespace Optima {
namespace {

struct WarningState
{
    std::mutex mutex;
    std::shared_ptr<const WarningHandler> handler; // null: use the default (print to std::cout)
};

auto warningState() -> WarningState&
{
    static WarningState state;
    return state;
}

} // namespace

auto setWarningHandler(WarningHandler handler) -> void
{
    // Build the new handler before locking, and let the displaced one be destroyed after unlocking: destroying it can
    // run user code (the destructor of something the handler captured) that calls setWarningHandler again, which
    // would deadlock on the mutex if it ran inside the lock.
    std::shared_ptr<const WarningHandler> replacement =
        handler ? std::make_shared<const WarningHandler>(std::move(handler)) : nullptr;

    auto& state = warningState();
    {
        std::lock_guard<std::mutex> lock(state.mutex);
        state.handler.swap(replacement);
    }
    // `replacement` now holds the previous handler; it is released here, outside the lock
}

namespace internal {

auto emitWarning(const std::string& message) -> void
{
    std::shared_ptr<const WarningHandler> handler;
    {
        auto& state = warningState();
        std::lock_guard<std::mutex> lock(state.mutex);
        handler = state.handler;
    }
    if(handler)
        (*handler)(message);   // outside the lock: the handler may itself log or call setWarningHandler
    else
        std::cout << "\033[1;33m***OPTIMA WARNING***\033[0m " << message << "\n";
}

} // namespace internal
} // namespace Optima
