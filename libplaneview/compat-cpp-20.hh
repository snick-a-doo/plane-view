// Copyright © 2022 Sam Varner
//
// This file is part of Plane View
//
// Plane View is free software: you can redistribute it and/or modify it under the terms
// of the GNU General Public License as published by the Free Software Foundation, either
// version 3 of the License, or (at your option) any later version.
//
// Plane View is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with Plane
// View.  If not, see <http://www.gnu.org/licenses/>.

/// Definitions introduced in C++ 20
#if __cplusplus >= 202000L
#error "To be used only when compiling to pre C++ 20 standards."
#endif
namespace compat
{
    /// Replacement for std::numbers::pi.
    namespace numbers
    {
        double constexpr pi{3.14159265358979323846};
    }

    /// Replacement for std::midpoint(). Average of a and b without overflow.
    template <typename T> T midpoint(T a, T b)
    {
        return a + (b - a)/2;
    }
}
