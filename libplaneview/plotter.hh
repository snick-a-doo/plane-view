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
    Line_Style m_line_style{Line_Style::points};
};

#endif // PLANE_VIEW_LIBPLANEVIEW_GRID_MAP_HH_INCLUDED
