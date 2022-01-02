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
#include <deque>
#include <string>
#include <vector>

/// The Cairo drawing context.
using Context = Cairo::RefPtr<Cairo::Context>;
/// Alias for a vector of vectors.
using VV = std::vector<V>;

enum class Line_Style
{
    points = 0,
    lines,
    lines_and_points,
    num,
};

/// A drawing area that show a graph and handles interaction.
class Plotter : public Gtk::DrawingArea
{
public:
    Plotter(Glib::RefPtr<Gtk::Application> app);

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

    /// Callback for reading standard input.
    bool on_read(Glib::IOCondition io_cond);

    /// Set the axes wide enough to contain all points plus padding, then redraw.
    void autoscale();

    /// A vector of x-value vectors, one for each trace.
    VV m_xss;
    /// A vector of y-value vectors, one for each trace.
    VV m_yss;
    /// The style of points and lines on the graph.
    Line_Style m_line_style{Line_Style::points};
    /// The horizontal axis object.
    Axis m_x_axis;
    /// The vertical axis object.
    Axis m_y_axis;

    /// The pixels where a dragging started. No value if dragging is not in progress.
    std::optional<double> m_drag_start_x{0.0};
    std::optional<double> m_drag_start_y{0.0};
    /// The pointer position while dragging.
    double m_drag_x;
    double m_drag_y;

    /// The input channel for reading data.
    Glib::RefPtr<Glib::IOChannel> m_io_channel;
    /// The signal handler for reading data. Disconnected when reading finishes.
    sigc::connection m_io_connection;

    /// The application object. Used for calling quit().
    Glib::RefPtr<Gtk::Application> m_app;

    /// Undo
    /// @{
    /// State information handled by the undo stack.
    struct State
    {
        double x_min; ///< X-axis low range.
        double x_max; ///< X-axis high range.
        double y_min; ///< Y-axis low range.
        double y_max; ///< Y-axis high range.
    };
    /// Add a state to the undo stack.
    void record();
    /// Move to the previous state in the undo stack.
    void undo();
    /// Move to the next state in the undo stack.
    void redo();
    /// Update the view to match the current position in the undo stack.
    void update(std::deque<State>::const_iterator it);
    /// The undo stack.
    std::deque<State> m_history;
    /// The current state in the history.
    std::deque<State>::const_iterator m_now;
    /// @}
};

#endif // PLANE_VIEW_LIBPLANEVIEW_GRID_MAP_HH_INCLUDED
