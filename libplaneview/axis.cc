// Copyright © 2021-2022 Sam Varner
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

#include "axis.hh"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <stdexcept>

const int min_ticks{4};

// Given the range of an axis and the minimum number of tick marks, return a rounded
// interval between ticks and the maximum number of non-zero digits after the decimal for
// multiples of that range.
static std::pair<double, int> axis_round(double range, double ticks)
{
    using namespace std;

    // The unrounded interval.
    auto interval{std::abs(range)/ticks};
    // The place of the most significant digit of the interval.
    auto exponent{static_cast<int>(floor(log10(interval)))};
    // The interval shifted to the range [1, 10)
    auto mantissa{interval*pow(10, -exponent)};
    auto closest{100.0};
    auto last_m{100.0};
    auto prec{1};
    for (auto m : {1.0, 2.0, 2.5, 5.0, 10.0})
    {
        auto diff{abs(mantissa - m)};
        if (diff >= closest)
            break;
        closest = diff;
        last_m = m;
    }
    if (last_m == 100)
        ++exponent;

    auto rounded{last_m*pow(10, exponent)};
    prec = (last_m == 2.5 ? 1 : last_m == 10.0 ? -1 : 0) - exponent;
    return {rounded, prec};
}

void Axis::set_pos(double low_pos, double high_pos)
{
    m_low_pos = low_pos;
    m_high_pos = high_pos;
}

void Axis::set_coord_range(double low_coord, double high_coord)
{
    auto cant_be_int = [](double x) {
        return std::isnan(x) || std::isinf(x)
            || x < std::numeric_limits<int>::min()
            || x > std::numeric_limits<int>::max();
    };
    auto [dx, prec] = axis_round(high_coord - low_coord, min_ticks);
    // Double the tick density to include minor ticks.
    dx *= 0.5;
    // See if we can handle the range being asked for.
    if (std::abs(dx) < std::numeric_limits<double>::epsilon()
        || cant_be_int(low_coord/dx) || cant_be_int(high_coord/dx))
    {
        std::cerr << "Range out of range: "
                  << low_coord << ' ' << high_coord << " dx = " << dx << std::endl;
        return;
    }
    m_low_coord = low_coord;
    m_high_coord = high_coord;
}

void Axis::set_coord_range_by_pos(double low_pos, double high_pos)
{
    auto coord1{pos_to_coord(low_pos)};
    auto coord2{pos_to_coord(high_pos)};
    set_coord_range(std::min(coord1, coord2), std::max(coord1, coord2));
}

void Axis::move_coord_range_by_pos(double delta_pos)
{
    auto delta_coord{delta_pos*(m_high_coord - m_low_coord)/(m_high_pos - m_low_pos)};
    set_coord_range(m_low_coord + delta_coord, m_high_coord + delta_coord);
}

void Axis::scale_range(double factor, std::optional<double> center_pos)
{
    auto mid{center_pos
        ? pos_to_coord(*center_pos)
        : std::midpoint(m_low_coord, m_high_coord)};
    set_coord_range(factor*(m_low_coord - mid) + mid, factor*(m_high_coord - mid) + mid);
}

bool is_in_range(double x, double a, double b)
{
    // Compare in a way doesn't care which endpoint is higher.
    return (a < x) == (x < b) || x == a || x == b;
}

bool Axis::is_in_pos_range(double pos) const
{
    return is_in_range(pos, m_low_pos, m_high_pos);
}

bool Axis::is_in_coord_range(double coord) const
{
    return is_in_range(coord, m_low_coord, m_high_coord);
}

std::pair<double, double> Axis::get_coord_range() const
{
    return {m_low_coord, m_high_coord};
}

std::pair<double, double> Axis::get_pos_range() const
{
    return {m_low_pos, m_high_pos};
}

double Axis::pos_to_coord(double pos) const
{
    auto scale{(pos - m_low_pos)/(m_high_pos - m_low_pos)};
    return std::lerp(m_low_coord, m_high_coord, scale);
}

double Axis::coord_to_pos(double coord) const
{
    auto scale{(coord - m_low_coord)/(m_high_coord - m_low_coord)};
    return std::lerp(m_low_pos, m_high_pos, scale);
}

std::vector<double> Axis::coord_to_pos(std::vector<double> const& coords) const
{
    std::vector<double> poss;
    for (auto coord : coords)
        poss.push_back(coord_to_pos(coord));
    return poss;
}

std::vector<Axis::Tick> Axis::get_ticks() const
{
    std::vector<Axis::Tick> ts;
    auto const dx{0.5 * axis_round(m_high_coord - m_low_coord, min_ticks).first};
    // Double the tick density to include minor ticks.
    auto low{static_cast<int>(std::ceil(m_low_coord/dx))};
    auto high{static_cast<int>(std::floor(m_high_coord/dx))};

    for (auto x{low}; x <= high; ++x)
    {
        auto coord{x*dx};
        ts.emplace_back(coord_to_pos(coord), std::nullopt);
        // Even-numbered ticks are major.
        if (x % 2 == 0)
            ts.back().label = format(coord);
    }
    return ts;
}

double Axis::round_pos(double pos, int divisions) const
{
    auto dx{axis_round(m_high_coord - m_low_coord, divisions).first};
    auto coord{pos_to_coord(pos)};
    return coord_to_pos(dx * std::round(coord/dx));
}

std::string Axis::format(double x, int extra_prec) const
{
    auto prec{axis_round(m_high_coord - m_low_coord, min_ticks).second};
    std::ostringstream os;
    os << std::fixed << std::setprecision(std::max(0, prec + extra_prec)) << x;
    return os.str();
}
