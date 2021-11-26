// Copyright © 2021 Sam Varner
//
// This file is part of Plane View.
//
// 4color is free software: you can redistribute it and/or modify it under the terms of
// the GNU General Public License as published by the Free Software Foundation, either
// version 3 of the License, or (at your option) any later version.
//
// 4color is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with 4color.
// If not, see <http://www.gnu.org/licenses/>.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <plotter.hh>

TEST_CASE("ticks")
{
    Axis ax;
    ax.set_pixels(100, 1000);
    ax.set_range(1.0, 10.0);
    auto ts{ax.ticks()};
    CHECK(ts.size() == 5);
    CHECK(ts[0].pixel == 200);
    CHECK(ts[0].label == "2");
    CHECK(ts[1].pixel == 400);
    CHECK(ts[1].label == "4");
    CHECK(ts[4].pixel == 1000);
    CHECK(ts[4].label == "10");
}
