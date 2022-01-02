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

#ifndef PLANE_VIEW_LIBPLANEVIEW_AXIS_HH_INCLUDED
#define PLANE_VIEW_LIBPLANEVIEW_AXIS_HH_INCLUDED

#include <string>
#include <vector>

using V = std::vector<double>;

class Axis
{
public:
    Axis() = default;
    ~Axis() = default;

    void set_pixels(int low, int high);
    void set_label_pos(int pos);
    void set_range(double low, double high);
    void set_range_pixels(double low, double high);

    /// @return A pair with the endpoints of the axis.
    std::pair<double, double> get_range() const;
    int low_pos() const { return m_low_pos; }
    int high_pos() const { return m_high_pos; }
    int label_pos() const { return m_label_pos; }
    int size() const { return m_high_pos - m_low_pos; }
    V to_pixels(V const& xs) const;
    double to_coord(double x) const;
    void zoom(double factor);
    std::string format(double x, int extra_prec = 0) const;

    struct Point
    {
        double pixel;
        std::string label;
    };
    using VPoint = std::vector<Point>;
    VPoint ticks() const;

private:
    double to_pixels(double x) const;
    int m_low_pos{0};
    int m_high_pos{100};
    int m_label_pos{0};
    double m_min{0.0};
    double m_max{1.0};
    mutable int m_precision{0};
};

#endif // PLANE_VIEW_LIBPLANEVIEW_AXIS_HH_INCLUDED
