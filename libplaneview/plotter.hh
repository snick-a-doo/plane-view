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

#ifndef PLANE_VIEW_LIBPLANEVIEW_GRID_MAP_HH_INCLUDED
#define PLANE_VIEW_LIBPLANEVIEW_GRID_MAP_HH_INCLUDED

#include "axis.hh"

#include <gtkmm.h>

#include <optional>
#include <string>
#include <vector>

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
    Plotter(Glib::RefPtr<Gtk::Application> app);
    ~Plotter();
    void plot(V const& xs, V const& ys);

private:
    /// DrawingArea methods
    /// @{
    virtual bool on_key_press_event(GdkEventKey* event) override;
    virtual bool on_button_press_event(GdkEventButton* event) override;
    virtual bool on_motion_notify_event(GdkEventMotion* event) override;
    virtual bool on_button_release_event(GdkEventButton* event) override;
    /// Callback for size change.
    virtual bool on_configure_event(GdkEventConfigure* event) override;
    virtual bool on_draw(Context const& cr) override;
    /// @}

    bool on_read(Glib::IOCondition io_cond);

    void autoscale();

    V m_xs;
    V m_ys;
    Line_Style m_line_style{Line_Style::points};
    Axis m_x_axis;
    Axis m_y_axis;

    std::string m_pipe;
    Glib::RefPtr<Glib::IOChannel> m_io_channel;
    V m_read_xs;
    V m_read_ys;
    int m_read_state{-1};

    Glib::RefPtr<Gtk::Application> m_app;

    std::optional<double> m_drag_start_x{0.0};
    std::optional<double> m_drag_start_y{0.0};
    double m_drag_x;
    double m_drag_y;
};

#endif // PLANE_VIEW_LIBPLANEVIEW_GRID_MAP_HH_INCLUDED
