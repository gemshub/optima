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

#pragma once

// C++ includes
#include <exception>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>

namespace Optima {
namespace internal {

template <typename Arg>
auto _stringfy(std::stringstream& ss, const std::string& sep, const Arg& item) -> void
{
    ss << item;
}

template <typename Arg, typename... Args>
auto _stringfy(std::stringstream& ss, const std::string& sep, const Arg& item, Args... items) -> void
{
    ss << item << sep;
    _stringfy(ss, sep, items...);
}

/// Concatenate the arguments into a string using a given separator string.
template <typename... Args>
auto stringfy(const std::string& sep, Args... items) -> std::string
{
    std::stringstream ss;
    _stringfy(ss, sep, items...);
    return ss.str();
}

template <typename... Args>
auto str(Args... items) -> std::string
{
    return stringfy("", items...);
}

} // namespace internal

/// The type of a function that receives the text of an Optima warning (without the "OPTIMA WARNING" tag).
using WarningHandler = std::function<void(const std::string&)>;

/// Replace what happens to Optima warnings. The default handler prints them to `std::cout`.
/// Pass an empty function to restore the default; pass `[](const std::string&){}` to silence them.
/// The handler may be called from any thread that runs the solver.
auto setWarningHandler(WarningHandler handler) -> void;

namespace internal {

/// Deliver a warning message to the current warning handler.
auto emitWarning(const std::string& message) -> void;

} // namespace internal

/// Issue a warning message if condition is true.
template<typename... Args>
auto warning(bool condition, Args... items) -> void
{
    if(condition)
        internal::emitWarning(internal::str(items...));
}

/// Raise a runtime error if condition is true.
template<typename... Args>
auto error(bool condition, Args... items) -> void
{
    if(condition)
        throw std::runtime_error(internal::str("\033[1;31m***OPTIMA ERROR***\033[0m\n", internal::str(items...), "\n"));
}

/// Define a macro to print a warning messageif condition is true.
/// @warning Note the use of ... and __VA_ARGS__ in the implementation.
/// @warning The use of `#define macro(args...) function(args)` causes error in
/// @warning MSVC compilers if compiler option `/Zc:preprocessor` (available in Visual
/// @warning Studio 2019 v16.5) is not specified!
/// @ingroup Common
#define warningif(condition, ...) \
    { \
        if((condition)) { \
            ::Optima::internal::emitWarning(internal::str(__VA_ARGS__)); \
        } \
    }

/// Define a macro to raise a runtime exception if condition is true.
/// @warning Note the use of ... and __VA_ARGS__ in the implementation.
/// @warning The use of `#define macro(args...) function(args)` causes error in
/// @warning MSVC compilers if compiler option `/Zc:preprocessor` (available in Visual
/// @warning Studio 2019 v16.5) is not specified!
/// @ingroup Common
#define errorif(condition, ...) \
    { \
        if((condition)) { \
            throw std::runtime_error(internal::str("\033[1;31m***OPTIMA ERROR***\033[0m\n", internal::str(__VA_ARGS__), "\n")); \
        } \
    }

} // namespace Optima
