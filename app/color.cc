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

#include <color.hh>

#include <array>
#include <cassert>
#include <cmath>
#include <numbers>

/// The Brewer Set 1 color palette.
std::array<Color, 8> constexpr brewer_1_colors {
    Color{0xe4, 0x1a, 0x1c}, Color{0x37, 0x7e, 0xb8},
    Color{0x4d, 0xaf, 0x4a}, Color{0x98, 0x4e, 0xa3},
    Color{0xff, 0x7f, 0x00}, Color{0xff, 0xff, 0x33},
    Color{0xa6, 0x56, 0x28}, Color{0xf7, 0x81, 0xbf},
};

/// The Brewer Set 2 color palette.
std::array<Color, 8> constexpr brewer_2_colors {
    Color{0x66, 0xc2, 0xa5}, Color{0xfc, 0x8d, 0x62},
    Color{0x8d, 0xa0, 0xcb}, Color{0xe7, 0x8a, 0xc3},
    Color{0xa6, 0xd8, 0x54}, Color{0xff, 0xd9, 0x2f},
    Color{0xe5, 0xc4, 0x94}, Color{0xb3, 0xb3, 0xb3},
};

/// The Okabe-Ito color palette.
std::array<Color, 8> constexpr okabe_ito_colors{
    Color{230, 159, 0}, Color{86, 180, 223},
    Color{0, 158, 115}, Color{240, 228, 66},
    Color{0, 114, 178}, Color{213, 94, 0},
    Color{204, 121, 167}, Color{0, 0, 0}
};

/// Hue, chroma, lightness to RGB conversion from R source code.
/// https://github.com/wch/r-source/blob/trunk/src/library/grDevices/src/colors.c
/// @param h Hue in degrees from red.
/// @param c Chroma -- color intensity, 0 to 100.
/// @param v Value -- Lightness, 0 to 100.
/// @return A Color object with integer RGB from 0 to 255.
Color R_HCL_to_RGB(double h, double c, double l)
{
    if (l <= 0.0)
        return {};

    // Convert to CIE-LUV.
    h *= std::numbers::pi/180.0;
    auto L{l};
    auto U{c*cos(h)};
    auto V{c*sin(h)};

    // Convert to CIE-XYZ.
    auto X{0.0};
    auto Y{0.0};
    auto Z{0.0};
    if (L > 0.0 || U != 0.0 || V != 0.0)
    {
        Y = 100.0*(L > 7.999592 ? std::pow((L + 16.0)/116.0, 3) : L/903.3);
        auto u{U/(13.0*L) + 0.1978398};
        auto v{V/(13.0*L) + 0.4683363};
        X = 9.0*Y*u/(4.0*v);
        Z = -X/3.0 - 5.0*Y + 3.0*Y/v;
    }

    // Convert CIE-XYZ to sRGB.
    auto gtrans = [](double u) {
        auto constexpr gamma{2.4};
        return static_cast<int>(
            255*(u > 0.00304 ? 1.055*std::pow(u, 1.0/gamma) - 0.055 : 12.92*u));
    };

    return {gtrans(0.01*( 3.240479*X - 1.537150*Y - 0.498535*Z)),
            gtrans(0.01*(-0.969256*X + 1.875992*Y + 0.041556*Z)),
            gtrans(0.01*( 0.055648*X - 0.204043*Y + 1.057311*Z))};
}

Color get_color(Palette palette, std::size_t index, std::size_t total)
{
    assert(total != 0);
    switch (palette)
    {
    case Palette::ggplot2:
        // Equally spaced hues starting from 15 degrees. Constant chroma and lightness.
        return R_HCL_to_RGB(std::lerp(15, 375, double(index)/total), 100, 65);
    case Palette::brewer_set_1:
        return brewer_1_colors[index % brewer_1_colors.size()];
    case Palette::brewer_set_2:
        return brewer_2_colors[index % brewer_2_colors.size()];
    case Palette::okabe_ito:
        return okabe_ito_colors[index % okabe_ito_colors.size()];
    default:
        assert(false);
    }
}
