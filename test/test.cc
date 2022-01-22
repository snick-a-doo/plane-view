// Copyright © 2021 Sam Varner
//
// This file is part of Plane View.
//
// 4color is free software: you can redistribute it and/or modify it under the terms of
// the GNU General Public License as published by the Free Software Foundation, either
// version 3 of the License, or (at your option) any later version.
//
// Plane View is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with Plane
// View.  If not, see <http://www.gnu.org/licenses/>.

#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include "doctest.h"

#include <plotter.hh>

bool close(double x1, double x2, double tolerance = 1e-9)
{
    if (std::abs(x1 - x2) < tolerance)
        return true;
    std::cerr << x1 << " != " << x2 << ' ' << x1 - x2 << std::endl;
    return false;
}

TEST_CASE("ticks 10")
{
    Axis ax;
    ax.set_pos(0, 1000, 10);
    ax.set_coord_range(1.0, 4.0, 0.05);
    // [1,4] is padded to [0.85,4.15]. Not affected by rounding
    auto [low, high] = ax.get_coord_range();
    CHECK(low == 0.85);
    CHECK(high == 4.15);
    auto slope{1000.0/3.3};
    auto ts{ax.get_ticks()};
    CHECK(ts.size() == 4);
    CHECK(close(ts[0].position, 0.15*slope));
    CHECK(close(ts[1].position - ts[0].position, ts[2].position - ts[1].position));
    CHECK(close(ts[3].position - ts[2].position, ts[2].position - ts[1].position));
    CHECK(close(ts[3].position, 1000 - 0.15*slope));

    CHECK(ts[0].label == "1");
    CHECK(ts[1].label == "2");
    CHECK(ts[2].label == "3");
    CHECK(ts[3].label == "4");
}

TEST_CASE("ticks 10 round")
{
    Axis ax;
    ax.set_pos(0, 1000, 10);
    ax.set_coord_range(1.01234, 4.05678, 0.05);
    // [1,4] is padded to [0.860118,4.209002] and rounded to [0.86,4.21]
    auto [low, high] = ax.get_coord_range();
    CHECK(low == 0.86);
    CHECK(high == 4.21);
    auto slope{1000.0/3.35};
    auto ts{ax.get_ticks()};
    CHECK(ts.size() == 4);
    CHECK(close(ts[0].position, 0.14*slope));
    CHECK(close(ts[1].position - ts[0].position, ts[2].position - ts[1].position));
    CHECK(close(ts[3].position - ts[2].position, ts[2].position - ts[1].position));
    CHECK(close(ts[3].position, 1000 - 0.21*slope));

    CHECK(ts[0].label == "1");
    CHECK(ts[1].label == "2");
    CHECK(ts[2].label == "3");
    CHECK(ts[3].label == "4");
}

TEST_CASE("ticks 2")
{
    Axis ax;
    ax.set_pos(0, 1000, 10);
    ax.set_coord_range(1.0, 8.0, 0.05);
    auto ts{ax.get_ticks()};
    CHECK(ts.size() == 4);
    CHECK(ts[0].label == "2");
    CHECK(ts[1].label == "4");
    CHECK(ts[2].label == "6");
    CHECK(ts[3].label == "8");
}

TEST_CASE("ticks 2.5")
{
    Axis ax;
    ax.set_pos(0, 1000, 10);
    ax.set_coord_range(1.0, 10.0, 0.05);
    auto ts{ax.get_ticks()};
    CHECK(ts.size() == 4);
    CHECK(ts[0].label == "2.5");
    CHECK(ts[1].label == "5.0");
    CHECK(ts[2].label == "7.5");
    CHECK(ts[3].label == "10.0");
}

TEST_CASE("ticks 5")
{
    Axis ax;
    ax.set_pos(0, 1000, 10);
    ax.set_coord_range(1.0, 20.0, 0.05);
    // [1,20] is padded to [0.05,20.95]. Not affected by rounding.
    auto [low, high] = ax.get_coord_range();
    CHECK(low == 0.05);
    CHECK(high == 20.95);
    auto ts{ax.get_ticks()};
    CHECK(ts.size() == 4);
    CHECK(ts[0].label == "5");
    CHECK(ts[1].label == "10");
    CHECK(ts[2].label == "15");
    CHECK(ts[3].label == "20");
}

TEST_CASE("ticks 1")
{
    Axis ax;
    ax.set_pos(0, 1000, 10);
    ax.set_coord_range(1.0, 5.0, 0.05);
    auto ts{ax.get_ticks()};
    CHECK(ts.size() == 5);
    CHECK(ts[0].label == "1");
    CHECK(ts[1].label == "2");
    CHECK(ts[2].label == "3");
    CHECK(ts[3].label == "4");
    CHECK(ts[4].label == "5");
}

TEST_CASE("ticks small")
{
    Axis ax;
    ax.set_pos(0, 1000, 10);
    ax.set_coord_range(1.0, 1.01, 0.05);
    auto [low, high] = ax.get_coord_range();
    CHECK(low == 0.9995);
    CHECK(high == 1.0105);
    auto ts{ax.get_ticks()};
    CHECK(ts.size() == 5);
    CHECK(ts[0].label == "1.0000");
    CHECK(ts[1].label == "1.0025");
    CHECK(ts[2].label == "1.0050");
    CHECK(ts[3].label == "1.0075");
    CHECK(ts[4].label == "1.0100");
}

TEST_CASE("ticks big")
{
    Axis ax;
    ax.set_pos(0, 1000, 10);
    ax.set_coord_range(200.0, 2000, 0.05);
    auto ts{ax.get_ticks()};
    CHECK(ts.size() == 4);
    CHECK(ts[0].label == "500");
    CHECK(ts[1].label == "1000");
    CHECK(ts[2].label == "1500");
    CHECK(ts[3].label == "2000");
}

TEST_CASE("closest")
{
    Axis xax;
    xax.set_pos(100, 1000, 10);
    xax.set_coord_range(20, 40);
    Axis yax;
    yax.set_pos(200, 2, 10);
    yax.set_coord_range(-20, 80);

    SUBCASE("empty")
    {
        V xs;
        V ys;
        CHECK(!find_closest_point({200, 100}, xs, ys, xax, yax));
    }
    SUBCASE("one below range")
    {
        V xs{10};
        V ys{1};
        CHECK(!find_closest_point({200, 100}, xs, ys, xax, yax));
    }
    SUBCASE("one in range")
    {
        V xs{30};
        V ys{1};
        CHECK(!find_closest_point({0, 0}, xs, ys, xax, yax));
        CHECK(!find_closest_point({0, 100}, xs, ys, xax, yax));
        auto p{find_closest_point({200, 100}, xs, ys, xax, yax)};
        CHECK(p);
        CHECK(*p == Point(30, 1));
        p = find_closest_point({900, 100}, xs, ys, xax, yax);
        CHECK(p);
        CHECK(*p == Point(30, 1));
        CHECK(!find_closest_point({1001, 100}, xs, ys, xax, yax));
    }
    SUBCASE("two in range")
    {
        V xs{0, 15, 30, 45};
        V ys{20, 20, 20, 20};
        auto p{find_closest_point({101, 10}, xs, ys, xax, yax)};
        CHECK(p);
        CHECK(*p == Point(15, 20));
        p = find_closest_point({500, 10}, xs, ys, xax, yax);
        CHECK(p);
        CHECK(*p == Point(30, 20));
        p = find_closest_point({900, 10}, xs, ys, xax, yax);
        CHECK(p);
        CHECK(*p == Point(45, 20));
    }
}
