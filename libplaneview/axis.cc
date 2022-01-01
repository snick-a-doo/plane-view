// Copyright © 2021 Sam Varner
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
std::pair<double, int> axis_round(double range, double ticks)
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
    prec = std::max(0, (last_m == 2.5 ? 1 : last_m == 10.0 ? -1 : 0) - exponent);
    return {rounded, prec};
}

void Axis::set_pixels(int low, int high)
{
    m_low_pos = low;
    m_high_pos = high;
}

void Axis::set_range(double low, double high)
{
    auto pad{0.05*(high - low)};
    m_min = low - pad;
    m_max = high + pad;
}

void Axis::set_range_pixels(double low, double high)
{
    auto r1{to_coord(low)};
    auto r2{to_coord(high)};
    m_min = std::min(r1, r2);
    m_max = std::max(r1, r2);
}

void Axis::zoom(double factor)
{
    auto mid{std::midpoint(m_min, m_max)};
    auto dist{0.5*(m_max - m_min)/factor};
    m_min = mid - dist;
    m_max = mid + dist;
}

std::pair<double, double> Axis::get_range() const
{
    return {m_min, m_max};
}

void Axis::set_label_pos(int pos)
{
    m_label_pos = pos;
}

Axis::VPoint Axis::ticks() const
{
    VPoint ts;
    auto [dx, x_prec] = axis_round((m_max - m_min), min_ticks);
    auto low{static_cast<int>(std::ceil(m_min/dx))};
    auto high{static_cast<int>(std::floor(m_max/dx))};
    for (auto x{low}; x <= high; ++x)
    {
        std::ostringstream os;
        os << std::fixed << std::setprecision(x_prec) << x*dx;
        ts.emplace_back(to_pixels(x*dx), os.str());
    }
    return ts;
}

V Axis::to_pixels(V const& xs) const
{
    V px;
    for (auto x : xs)
        px.push_back(to_pixels(x));
    return px;
}

double Axis::to_pixels(double x) const
{
    auto scale{(m_high_pos - m_low_pos)/(m_max - m_min)};
    return m_low_pos + scale*(x - m_min);
}

double Axis::to_coord(double x) const
{
    auto scale{(m_high_pos - m_low_pos)/(m_max - m_min)};
    return (x - m_low_pos)/scale + m_min;
}
