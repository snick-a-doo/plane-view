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

#ifndef PLANE_VIEW_LIBPLANEVIEW_COLOR_HH_INCLUDED
#define PLANE_VIEW_LIBPLANEVIEW_COLOR_HH_INCLUDED

#include <cstddef>

/// RGBA color struct.
struct Color
{
    int red{0};
    int green{0};
    int blue{0};
    int alpha{255}; ///< Opaque by default.
};

/// The available color sets.
enum class Palette
{
    ggplot2,
    brewer_set_1,
    brewer_set_2,
    okabe_ito,
};

/// @param palette One of the available color sets.
/// @pacam index Which color in the set.
/// @param total The number of colors in the set, if not fixed. Currently only used with
/// the ggplot2 palette.
/// @return A color struct from the specified palette.
Color get_color(Palette palette, std::size_t index, std::size_t total);

#endif // PLANE_VIEW_LIBPLANEVIEW_COLOR_HH_INCLUDED
