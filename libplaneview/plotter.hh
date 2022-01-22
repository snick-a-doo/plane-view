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

#ifndef PLANE_VIEW_LIBPLANEVIEW_GRID_MAP_HH_INCLUDED
#define PLANE_VIEW_LIBPLANEVIEW_GRID_MAP_HH_INCLUDED

#include "axis.hh"

#include <gtkmm.h>

#include <optional>
#include <deque>
#include <string>
#include <vector>

/// A 2D point.
struct Point
{
    double x{0};
    double y{0};
    // Generate comparison operators.
    auto operator <=>(Point const& point) const = default;
};
std::ostream& operator<<(std::ostream& os, Point const& p);

/// The Cairo drawing context.
using Context = Cairo::RefPtr<Cairo::Context>;
/// Alias for a vector of doubles.
using V = std::vector<double>;
/// Alias for a vector of vectors of doubles.
using VV = std::vector<V>;

enum class Line_Style
{
    points = 0,
    lines,
    lines_and_points,
    num,
};

/// @return The plot coordinates of the data point that is plotted closest to the
/// device-coordinate point p.
std::optional<Point> find_closest_point(Point const& p, V const& xs, V const& ys,
                                        Axis const& x_axis, Axis const& y_axis);

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
    virtual bool on_scroll_event(GdkEventScroll* event) override;
    /// Callback for size change.
    virtual bool on_configure_event(GdkEventConfigure* event) override;
    virtual bool on_draw(Context const& cr) override;
    /// @}

    /// Callback for reading standard input.
    bool on_read(Glib::IOCondition io_cond);

    /// Set the axes wide enough to contain all points plus padding, then redraw.
    void autoscale();

    /// Move the graph by the given fraction of the visible range. In the overview, move
    /// the subrange.
    void move(double x_frac, double y_frac);
    /// Scale the graph by the given fraction. In the overview, scale the subrange.
    void scale(double x_frac, double y_frac, std::optional<Point> center);

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

    enum class Side
    {
        left = 0, right, top, bottom, none,
    };

    /// Information about the range box in overview mode.
    class Subrange
    {
    public:
        /// Create a subrange with the given corners.
        Subrange(Point p1, Point p2);
        /// Change the corner positions.
        void set(Point p1, Point p2);
        /// Find which of the sides should change when the pointer moves depending on
        /// whether the point is close to a corner or edge, or in the interior.
        void start(Point p);
        /// Move the active sides
        void move(Point dp, bool all = false);
        /// Scale the subrange by the given fraction.
        void scale(double x_frac, double y_frac, std::optional<Point> center);
        /// @return The upper-left corner of the range, or if a side is specified, that
        /// side's active region.
        Point get_p1(Side side = Side::none) const;
        /// @return The lower-right corner of the range, or if a side is specified, that
        /// side's active region.
        Point get_p2(Side side = Side::none) const;
        /// @return The height of the range.
        double height() const;
        /// @return The width of the range.
        double width() const;

    private:
        /// Access side flags with a Side enum.
        bool& get_side(Side side);
        /// Flags that tell which sides of the range box are being adjusted.
        std::array<bool, 4> m_sides{false, false, false, false};
        /// The upper-left corner of the range in device coordinates.
        Point m_p1;
        /// The lower-right corner of the range in device coordinates.
        Point m_p2;
    };
    /// The active subrange in overview mode or nullptr.
    std::unique_ptr<Subrange> mp_subrange;

    /// Information about a drag.
    struct Drag
    {
        Point start; ///< Device position where the drag started.
        Point pointer; ///< Current pointer device position.
        bool shift{false}; ///< True if the shift key is pressed.
    };
    /// The current drag operation. No value if no drag is in progress.
    std::optional<Drag> m_drag;
    /// The coordinates of the data point that's closest to pointer. No value if it hasn't
    /// been requested or it was discarded.
    std::optional<Point> m_closest_point;

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
        bool incremental; ///< True for incremental changes.
        // Generate comparison operators.
        auto operator <=>(State const& state) const = default;
    };
    /// Add a state to the undo stack.
    /// @param incremental If true, the state is an incremental change. If the previous
    /// state was incremental, it will be overwritten.
    void record(bool incremental);
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
