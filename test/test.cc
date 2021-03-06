// Copyright © 2021-2022 Sam Varner
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


using Range = std::pair<double, double>;
bool range(Range r, double low, double high)
{
    return close(r.first, low) && close (r.second, high);
}

TEST_CASE("axis")
{
    Axis ax;
    ax.set_pos(10, 110);
    ax.set_coord_range(-1.0, 3.0);
    SUBCASE("initial")
    {
        CHECK(!ax.is_in_pos_range(1));
        CHECK(ax.is_in_pos_range(10));
        CHECK(ax.is_in_pos_range(50));
        CHECK(ax.is_in_pos_range(110));
        CHECK(!ax.is_in_pos_range(111));
        CHECK(!ax.is_in_coord_range(-2.0));
        CHECK(ax.is_in_coord_range(-1.0));
        CHECK(ax.is_in_coord_range(0.0));
        CHECK(ax.is_in_coord_range(3.0));
        CHECK(!ax.is_in_coord_range(3.00001));
        CHECK(range(ax.get_pos_range(), 10, 110));
        CHECK(range(ax.get_coord_range(), -1.0, 3.0));
    }
    SUBCASE("move")
    {
        ax.move_coord_range_by_pos(25.0);
        CHECK(range(ax.get_pos_range(), 10.0, 110.0)); // unchanged
        CHECK(range(ax.get_coord_range(), 0.0, 4.0));
    }
    SUBCASE("scale center")
    {
        ax.scale_range(-2.0);
        CHECK(range(ax.get_pos_range(), 10.0, 110.0)); // unchanged
        CHECK(range(ax.get_coord_range(), 5.0, -3.0));
    }
    SUBCASE("scale 3/4")
    {
        ax.scale_range(2.0, 85);
        CHECK(range(ax.get_pos_range(), 10.0, 110.0)); // unchanged
        CHECK(range(ax.get_coord_range(), -4.0, 4.0));
    }
}

TEST_CASE("inverted axis")
{
    Axis ax;
    ax.set_pos(10, 110);
    ax.set_coord_range(3.0, -1.0);
    SUBCASE("initial")
    {
        CHECK(!ax.is_in_pos_range(1));
        CHECK(ax.is_in_pos_range(10));
        CHECK(ax.is_in_pos_range(50));
        CHECK(ax.is_in_pos_range(110));
        CHECK(!ax.is_in_pos_range(111));
        CHECK(!ax.is_in_coord_range(-2.0));
        CHECK(ax.is_in_coord_range(-1.0));
        CHECK(ax.is_in_coord_range(0.0));
        CHECK(ax.is_in_coord_range(3.0));
        CHECK(!ax.is_in_coord_range(3.00001));
        CHECK(range(ax.get_pos_range(), 10, 110));
        CHECK(range(ax.get_coord_range(), 3.0, -1.0));
    }
    SUBCASE("move")
    {
        ax.move_coord_range_by_pos(25.0);
        CHECK(range(ax.get_pos_range(), 10.0, 110.0)); // unchanged
        CHECK(range(ax.get_coord_range(), 2.0, -2.0));
    }
    SUBCASE("scale center")
    {
        ax.scale_range(-2.0);
        CHECK(range(ax.get_pos_range(), 10.0, 110.0)); // unchanged
        CHECK(range(ax.get_coord_range(), -3.0, 5.0));
    }
    SUBCASE("scale 3/4")
    {
        ax.scale_range(2.0, 85);
        CHECK(range(ax.get_pos_range(), 10.0, 110.0)); // unchanged
        CHECK(range(ax.get_coord_range(), 6.0, -2.0));
    }
}

TEST_CASE("ticks 1 to 4")
{
    Axis ax;
    ax.set_pos(0, 1000);
    ax.set_coord_range(1.0, 4.0);
    auto [low, high] = ax.get_coord_range();
    CHECK(low == 1.0);
    CHECK(high == 4.0);
    auto slope{1000.0/3.0};
    auto ts{ax.get_ticks()};
    // Edge case: major ticks at 1.0, 1.5, 2.0,...4.0. Maybe we can say a tick at the
    // limit is worth 1/2. In that case 1,2,3,4 would only be 3 ticks and we need at least
    // 4.
    CHECK(ts.size() == 13);
    CHECK(close(ts[0].position, 0));
    CHECK(close(ts[1].position, 0.25*slope));
    CHECK(close(ts[12].position, 1000));

    CHECK(ts[0].label == "1.0");
    CHECK(!ts[1].label);
    CHECK(ts[2].label == "1.5");
    CHECK(!ts[3].label);
    CHECK(ts[4].label == "2.0");
    CHECK(!ts[5].label);
    CHECK(ts[12].label == "4.0");
}

TEST_CASE("ticks 1-e to 4+e")
{
    Axis ax;
    ax.set_pos(0, 1000);
    ax.set_coord_range(0.9999, 4.0001);
    auto [low, high] = ax.get_coord_range();
    CHECK(low == 0.9999);
    CHECK(high == 4.0001);
    auto slope{1000.0/3.0002};
    auto ts{ax.get_ticks()};
    CHECK(ts.size() == 7);
    CHECK(close(ts[0].position, 0.0001*slope));
    CHECK(close(ts[1].position - ts[0].position, ts[2].position - ts[1].position));
    CHECK(close(ts[3].position - ts[2].position, ts[2].position - ts[1].position));
    CHECK(close(ts[6].position, 1000 - 0.0001*slope));

    CHECK(ts[0].label == "1");
    CHECK(ts[2].label == "2");
    CHECK(ts[4].label == "3");
    CHECK(ts[6].label == "4");
}

TEST_CASE("ticks 2")
{
    Axis ax;
    ax.set_pos(0, 1000);
    ax.set_coord_range(1.0, 8.0);
    auto ts{ax.get_ticks()};
    CHECK(ts.size() == 8);
    CHECK(!ts[0].label);
    CHECK(ts[1].label == "2");
    CHECK(ts[3].label == "4");
    CHECK(ts[5].label == "6");
    CHECK(ts[7].label == "8");
}

TEST_CASE("ticks 2.5")
{
    Axis ax;
    ax.set_pos(0, 1000);
    ax.set_coord_range(1.0, 10.2);
    auto ts{ax.get_ticks()};
    CHECK(ts.size() == 8);
    CHECK(ts[1].label == "2.5");
    CHECK(ts[3].label == "5.0");
    CHECK(ts[5].label == "7.5");
    CHECK(ts[7].label == "10.0");
}

TEST_CASE("ticks 5")
{
    Axis ax;
    ax.set_pos(0, 1000);
    ax.set_coord_range(1.0, 20.0);
    auto ts{ax.get_ticks()};
    CHECK(ts.size() == 8);
    CHECK(ts[1].label == "5");
    CHECK(ts[3].label == "10");
    CHECK(ts[5].label == "15");
    CHECK(ts[7].label == "20");
}

TEST_CASE("ticks 1")
{
    Axis ax;
    ax.set_pos(0, 1000);
    ax.set_coord_range(1.0, 5.0);
    auto ts{ax.get_ticks()};
    CHECK(ts.size() == 9);
    CHECK(ts[0].label == "1");
    CHECK(ts[2].label == "2");
    CHECK(ts[4].label == "3");
    CHECK(ts[6].label == "4");
    CHECK(ts[8].label == "5");
}

TEST_CASE("ticks small")
{
    Axis ax;
    ax.set_pos(0, 1000);
    ax.set_coord_range(1.0, 1.01);
    auto ts{ax.get_ticks()};
    CHECK(ts.size() == 9);
    CHECK(ts[0].label == "1.0000");
    CHECK(ts[2].label == "1.0025");
    CHECK(ts[4].label == "1.0050");
    CHECK(ts[6].label == "1.0075");
    CHECK(ts[8].label == "1.0100");
}

TEST_CASE("ticks big")
{
    Axis ax;
    ax.set_pos(0, 1000);
    ax.set_coord_range(200.0, 2000);
    auto ts{ax.get_ticks()};
    CHECK(ts.size() == 8);
    CHECK(ts[1].label == "500");
    CHECK(ts[3].label == "1000");
    CHECK(ts[5].label == "1500");
    CHECK(ts[7].label == "2000");
}

TEST_CASE("closest")
{
    Axis xax;
    xax.set_pos(100, 1000);
    xax.set_coord_range(20, 40);
    Axis yax;
    yax.set_pos(200, 2);
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
