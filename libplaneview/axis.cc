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

const int min_ticks{4};

// Given the range of an axis and the minimum number of tick marks, return a rounded
// interval between ticks and the maximum number of non-zero digits after the decimal for
// multiples of that range.
static std::pair<double, int> axis_round(double range, double ticks)
{
    using namespace std;

    // The unrounded interval.
    auto interval{range/ticks};
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
    prec = (last_m == 2.5 ? 1 : last_m == 10.0 ? -1 : 0) - exponent;;
    return {rounded, prec};
}

/// Round the limits to one more digit of precision than the numbers on the tick marks.
static std::pair<double, double> rounded_range(double low_coord, double high_coord)
{
    // Set the precision for range limits to 2 more than the tick labels.
    auto [dx, prec] = axis_round(high_coord - low_coord, min_ticks);
    auto scale{std::pow(10, prec + 2)};
    auto round_prec = [scale](double x) {
        return std::round(x * scale) / scale;
    };
    return {round_prec(low_coord), round_prec(high_coord)};
}

void Axis::set_pos(double low_pos, double high_pos, double tick_label_pos)
{
    m_low_pos = low_pos;
    m_high_pos = high_pos;
    m_tick_label_pos = tick_label_pos;
}

void Axis::set_coord_range(double low_coord, double high_coord, double pad_fraction)
{
    auto extra{pad_fraction*(high_coord - low_coord)};
    auto [low, high] = rounded_range(low_coord - extra, high_coord + extra);
    m_low_coord = low;
    m_high_coord = high;
}

void Axis::set_pos_range(double low_pos, double high_pos, double pad_fraction)
{
    auto coord1{pos_to_coord(low_pos)};
    auto coord2{pos_to_coord(high_pos)};
    set_coord_range(std::min(coord1, coord2), std::max(coord1, coord2), pad_fraction);
}

void Axis::move_pos_range(double delta)
{
    auto scale{(m_high_pos - m_low_pos)/(m_high_coord - m_low_coord)};
    set_coord_range(m_low_coord + delta/scale, m_high_coord + delta/scale);
}

void Axis::scale_range(double factor, std::optional<double> center_pos)
{
    auto mid{center_pos
        ? pos_to_coord(*center_pos)
        : std::midpoint(m_low_coord, m_high_coord)};
    set_coord_range(factor*(m_low_coord - mid) + mid, factor*(m_high_coord - mid) + mid);
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
    auto scale{(m_high_pos - m_low_pos)/(m_high_coord - m_low_coord)};
    return (pos - m_low_pos)/scale + m_low_coord;
}

double Axis::coord_to_pos(double coord) const
{
    auto scale{(m_high_pos - m_low_pos)/(m_high_coord - m_low_coord)};
    return m_low_pos + scale*(coord - m_low_coord);
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
    auto [dx, prec] = axis_round((m_high_coord - m_low_coord), min_ticks);
    auto low{static_cast<int>(std::ceil(m_low_coord/dx))};
    auto high{static_cast<int>(std::floor(m_high_coord/dx))};
    for (auto x{low}; x <= high; ++x)
    {
        std::ostringstream os;
        os << std::fixed << std::setprecision(std::max(0, prec)) << x*dx;
        ts.emplace_back(coord_to_pos(x*dx), format(x*dx));
    }
    return ts;
}

std::string Axis::format(double x, int extra_prec) const
{
    auto [dx, prec] = axis_round(m_high_coord - m_low_coord, min_ticks);
    std::ostringstream os;
    os << std::fixed << std::setprecision(std::max(0, prec + extra_prec)) << x;
    return os.str();
}
