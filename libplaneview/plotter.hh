// Copyright © 2021 Sam Varner
//
// This file is part of Plane View
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

#ifndef PLANE_VIEW_LIBPLANEVIEW_GRID_MAP_HH_INCLUDED
#define PLANE_VIEW_LIBPLANEVIEW_GRID_MAP_HH_INCLUDED

#include <gtkmm.h>

#include <string>
#include <vector>

using V = std::vector<double>;
/// The Cairo drawing context.
using Context = Cairo::RefPtr<Cairo::Context>;

enum class Line_Style
{
    points = 0,
    lines,
    lines_and_points,
    num,
};

class Axis
{
public:
    Axis();
    void set_pixels(int low, int high);
    void set_label_pos(int pos);
    void set_range(double low, double high);
    int low_pos() const { return m_low_pos; }
    int high_pos() const { return m_high_pos; }
    int label_pos() const { return m_label_pos; }
    int size() const { return m_high_pos - m_low_pos; }
    V to_pixels(V const& xs) const;
    void zoom(double factor);

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
    const int m_min_ticks{5};
};

class Plotter : public Gtk::DrawingArea
{
public:
    Plotter();
    void plot(V const& xs, V const& ys);

private:
    /// DrawingArea methods
    /// @{
    virtual bool on_key_press_event(GdkEventKey* event) override;
    virtual bool on_button_press_event(GdkEventButton* event) override;
    /// Callback for size change.
    virtual bool on_configure_event(GdkEventConfigure* event) override;
    virtual bool on_draw(Context const& cr) override;
    /// @}

    V m_xs;
    V m_ys;
    double m_zoom{1.0};
    Line_Style m_line_style{Line_Style::points};
    Axis m_x_axis;
    Axis m_y_axis;
};

#endif // PLANE_VIEW_LIBPLANEVIEW_GRID_MAP_HH_INCLUDED
