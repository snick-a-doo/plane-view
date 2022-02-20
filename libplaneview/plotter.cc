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

#include <cassert>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <numbers>
#include <numeric>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <color.hh>
#include <plotter.hh>

/// The width of the border around the grid opposite the axes.
auto constexpr border_width{6.0};
/// The length of the tick marks outside of the grid.
auto constexpr tick_length{4.0};
auto constexpr text_size{10.0};
/// The size of the coordinates of the point near the pointer.
auto constexpr closest_point_text_size{18.0};
/// The length of the segments at the ends of the range bars.
auto constexpr range_bar_width{20.0};
/// The distance between the grid and the closest point on the numeric tick labels.
auto constexpr tick_label_gap{8.0};

/// The distance between the bottom of the grid and the bottom of the window.
auto constexpr grid_bottom_margin{border_width + closest_point_text_size + range_bar_width
    + text_size + tick_label_gap};
/// The distance between the left side of the grid and the left side of the window not
/// including the width of labels since that depends on the range.
auto constexpr grid_left_margin{border_width + range_bar_width + tick_label_gap};

/// The size of a plotted point.
auto constexpr point_radius{3.0};
/// The width of the adjustable border of the range box in overview mode.
auto constexpr handle_width{18.0};

/// The magnification for fine and coarse zooming.
std::pair constexpr zoom_factor{1.1, 2.0};
/// The fraction of the range to move for fine and coarse panning.
std::pair constexpr pan_distance{0.1, 1.0};
/// Number of divisions between 1st and last major tick. Cursor position is rounded to a
/// division boundary for coarse motion.
auto constexpr coarse_motion_divisions{20};
/// The padding as fraction of the range to add to an autoscaled graph.
auto constexpr autoscale_padding{0.05};

/// Start with a line style of line rather than points if there are more than this many
/// points. Also draw square points instead of round points.
std::size_t constexpr big_plot_threshold{5000};

/// The input string that indicates the start of a vector.
auto constexpr data_tag{"pv.data"};
/// The input string that indicates the end of data.
auto constexpr end_tag{"pv.end"};
/// The precision of the reported ranges relative to tick label precision.
auto constexpr output_precision{2};

/// An enumeration for the x and y directions.
enum class Direction {x, y};

///@{
/// Pre-defined colors.
Color constexpr black{0, 0, 0};
Color constexpr background{230, 230, 250}; // Tinted to avoid confusion with ggplot2.
Color constexpr white{255, 255, 255};
Color constexpr gray{128, 128, 128};
Color constexpr zoom_box_color{128, 128, 128, 64};
///}

enum class Point_Style {round, square};

std::ostream& operator<<(std::ostream& os, Point const& p)
{
    return os << p.x << ", " << p.y;
}

/// @return The argument restricted to the given range.
double clip(double x, double low, double high)
{
    return std::min(std::max(x, low), high);
};

/// Set Cairo's current active color using the Color struct above.
/// @param cr The Cairo graphics context.
/// @param color The RGBA color.
static void set_color(Context cr, Color const& color)
{
    // Cairo takes floating-point colors.
    cr->set_source_rgba(color.red/255.0,
                        color.green/255.0,
                        color.blue/255.0,
                        color.alpha/255.0);
}

/// Draw a rectangle given two points and a color.
/// @param p1 The upper-left corner.
/// @param p2 The lower-right corner.
/// @param color The fill or stroke color.
/// @param fill If true, the rectangle is filled, otherwise it's stroked.
static void draw_rectangle(Context cr, Point const& p1, Point const& p2,
                           Color const& color, bool fill)
{
    set_color(cr, color);
    cr->rectangle(p1.x, p1.y, p2.x - p1.x, p2.y - p1.y);
    if (fill)
        cr->fill();
    else
        cr->stroke();
}

/// Draw the axis numbers, tick marks, and grid lines.
static void draw_axes_and_grid(Context cr,
                               std::vector<Axis::Tick> const& x_ticks,
                               std::vector<Axis::Tick> const& y_ticks,
                               double x_low, double x_high,
                               double y_low, double y_high)
{
    // Draw the numbers and external tick lines.
    cr->set_line_width(1.0);
    set_color(cr, black);
    for (auto x : x_ticks)
    {
        if (!x.label)
            continue;
        Cairo::TextExtents ex;
        cr->get_text_extents(*x.label, ex);
        cr->move_to(x.position - 0.5*ex.width, y_low + tick_label_gap + ex.height);
        cr->show_text(*x.label);
        cr->move_to(x.position, y_low + tick_length);
        cr->line_to(x.position, y_low);
        cr->stroke();
    }
    for (auto y : y_ticks)
    {
        if (!y.label)
            continue;
        Cairo::TextExtents ex;
        cr->get_text_extents(*y.label, ex);
        cr->move_to(x_low - tick_label_gap - ex.width, y.position + 0.5*ex.height);
        cr->show_text(*y.label);
        cr->move_to(x_low - tick_length, y.position);
        cr->line_to(x_low, y.position);
        cr->stroke();
    }

    // Draw the grid.
    set_color(cr, background);
    cr->rectangle(x_low, y_low, x_high - x_low, y_high - y_low);
    cr->fill();
    set_color(cr, white);
    for (auto x : x_ticks)
    {
        cr->set_line_width(x.label ? 1.0 : 0.5);
        cr->move_to(x.position, y_low);
        cr->line_to(x.position, y_high);
        cr->stroke();
    }
    for (auto y : y_ticks)
    {
        cr->set_line_width(y.label ? 1.0 : 0.5);
        cr->move_to(x_low, y.position);
        cr->line_to(x_high, y.position);
        cr->stroke();
    }
}

/// Draw the data points and lines.
static void draw_plot(Context cr,
                      std::vector<double> const& pxs, std::vector<double> const& pys,
                      Color color, Point_Style point_style, Line_Style line_style,
                      std::optional<Point> point)
{
    auto N{std::min(pxs.size(), pys.size())};
    if (N == 0)
        return;

    set_color(cr, color);
    if (line_style != Line_Style::points)
    {
        cr->set_line_width(1.0);
        cr->move_to(pxs[0], pys[0]);
        for (std::size_t i = 1; i < N; ++i)
            cr->line_to(pxs[i], pys[i]);
        cr->stroke();
    }
    if (line_style != Line_Style::lines)
    {
        cr->set_line_width(2.0*point_radius);
        cr->set_line_cap(point_style == Point_Style::round
                         ? Cairo::LINE_CAP_ROUND : Cairo::LINE_CAP_SQUARE);
        for (std::size_t i = 0; i < N; ++i)
        {
            // Cairo doesn't draw square endcaps for zero-length lines because the
            // orientation is indeterminate. So draw a centipixel-long vertical line. For
            // a millipixel line, only some points are drawn.
            cr->move_to(pxs[i], pys[i]);
            cr->line_to(pxs[i], pys[i] + 0.01);
        }
        cr->stroke();
    }
    if (point)
    {
        // Draw crosshairs on the passed-in point.
        cr->set_line_width(0.5);
        set_color(cr, black);
        cr->move_to(point->x - 4*point_radius, point->y);
        cr->line_to(point->x + 4*point_radius, point->y);
        cr->move_to(point->x, point->y - 4*point_radius);
        cr->line_to(point->x, point->y + 4*point_radius);
        cr->stroke();
    }
}

/// Draw a range bar.
void draw_range(Context cr, int x0, int y0, int x1, int y1,
                std::string const& label, Direction axis)
{
    // Draw the range lines.
    auto dx{axis == Direction::x ? 0 : range_bar_width/2};
    auto dy{axis == Direction::x ? range_bar_width/2 : 0};
    set_color(cr, gray);
    cr->set_line_width(1.0);
    cr->move_to(x0 - dx, y0 - dy);
    cr->line_to(x0 + dx, y0 + dy);
    cr->move_to(x1 - dx, y1 - dy);
    cr->line_to(x1 + dx, y1 + dy);
    cr->move_to(x0, y0);
    cr->line_to(x1, y1);
    cr->stroke();

    // Draw the label.
    Cairo::TextExtents label_ext;
    cr->get_text_extents(label, label_ext);
    auto label_x{axis == Direction::x
        ? std::midpoint(x0, x1) - label_ext.width/2
        : border_width};
    auto label_y{axis == Direction::x
        ? y0 + label_ext.height/2
        : std::midpoint(y0, y1) + label_ext.height/2};
    auto gap{text_size - label_ext.height};
    cr->rectangle(label_x - gap, label_y - label_ext.height - gap,
                  label_ext.width + 2*gap, label_ext.height + 2*gap);
    set_color(cr, white);
    cr->fill();
    set_color(cr, black);
    cr->move_to(label_x, label_y);
    cr->show_text(label);
}

/// Draw range bars in the margins that correspond to the rectangle formed by p1 and p2.
void draw_ranges(Context cr, Glib::RefPtr<Gdk::Window> window,
                 Point const& p1, Point const& p2, Axis const& x_axis, Axis const& y_axis)
{
    auto dx{std::abs(x_axis.pos_to_coord(p2.x) - x_axis.pos_to_coord(p1.x))};
    auto label{x_axis.format(dx, 1)};
    auto y{window->get_height() - border_width - closest_point_text_size - range_bar_width/2};
    draw_range(cr, p1.x, y, p2.x, y, label, Direction::x);
    window->invalidate_rect(
        Gdk::Rectangle(0, y - range_bar_width/2, window->get_width(), range_bar_width), false);

    auto dy{std::abs(y_axis.pos_to_coord(p2.y) - y_axis.pos_to_coord(p1.y))};
    label = y_axis.format(dy, 1);
    auto x{border_width + range_bar_width/2};
    draw_range(cr, x, p1.y, x, p2.y, label, Direction::y);
    Cairo::TextExtents label_ext;
    cr->get_text_extents(label, label_ext);
    window->invalidate_rect(
        Gdk::Rectangle(border_width, 0,
                       std::max(label_ext.width + border_width, range_bar_width),
                       window->get_height()),
        false);

    // Show the slope.
    std::ostringstream os;
    os << "m = " << dy/dx;
    cr->move_to(border_width, y);
    cr->show_text(os.str());
}

/// Mark a an L-shaped region of the screen dirty after moving a corner of a rectangle.
/// @param p0 The position of the stationary corner, opposite the one that moved.
/// @param p1 The previous position of the corner that moved.
/// @param p2 The new position of the corner that moved.
/// @param width An extra amount added to each side of the region.
void redraw(Glib::RefPtr<Gdk::Window> window, Point p0, Point p1, Point p2, double width = 0)
{
    auto min3 = [](int a, int b, int c) {
        return std::min(a, std::min(b, c));
    };
    auto max3 = [](int a, int b, int c) {
        return std::max(a, std::max(b, c));
    };
    auto do_redraw = [window](double x1, double y1, double x2, double y2) {
        window->invalidate_rect(Gdk::Rectangle(x1, y1, x2-x1, y2-y1), false);
    };

    // Add 2 pixels to the width to avoid fringes.
    width += 2;

    // When a corner of the zoom box is moved, we need to redraw an L-shaped region
    // that has just been added or removed from the zoom box. First, mark the vertical
    // part of the L...
    do_redraw(std::min(p1.x, p2.x) - width, min3(p0.y, p1.y, p2.y) - width,
              std::max(p1.x, p2.x) + width, max3(p0.y, p1.y, p2.y) + width);
    // ...then the horizontal part.
    do_redraw(min3(p0.x, p1.x, p2.x) - width, std::min(p1.y, p2.y) - width,
              max3(p0.x, p1.x, p2.x) + width, std::max(p1.y, p2.y) + width);
}

/// @param p A point in device coordinates.
/// @param xs A vector of data coordinate x-values.
/// @param ys A vector of data coordinate y-values.
/// @param x_axis The x-axis abject.
/// @param y_axis The y-axis abject.
/// @return A point formed from the member of xs and the corresponding member of ys that's
/// plotted closest to p.
std::optional<Point> find_closest_point(Point const& p, V const& xs, V const& ys,
                                        Axis const& x_axis, Axis const& y_axis)
{
    // Ignore the point if it's outside the axes.
    if (!x_axis.is_in_pos_range(p.x) || !y_axis.is_in_pos_range(p.y))
        return std::nullopt;
    if (xs.empty() || ys.empty())
        return std::nullopt;

    // We want the point that appears closest to the pointer on the screen. Measure
    // distances in device coordinates.

    // Start with the 1st point whose x-position is not less than p.x. Use data
    // coordinates so we only have to convert the pointer position. That's okay because
    // we're only looking at the x-dimension.
    auto px{x_axis.pos_to_coord(p.x)};
    // Try 1st x that's greater than or equal to the point.
    auto it = std::find_if(xs.begin(), xs.end(), [px](double x) { return x >= px; });
    // If there's no such x, or it's off the graph, try the previous.
    if ((it == xs.end() || !x_axis.is_in_coord_range(*it)) && it != xs.begin())
        --it;
    // Give up if there's still no candidate.
    if (it == xs.end() || !x_axis.is_in_coord_range(*it))
        return std::nullopt;

    auto distance = [&](double dx, double y) {
        auto dy{y_axis.coord_to_pos(y) - p.y};
        return std::sqrt(dx*dx + dy*dy);
    };

    // Go through the x-values in order of closeness to the pointer in the
    // x-direction. Udate min_index when a closer point is found.
    std::size_t min_index{static_cast<std::size_t>(std::distance(xs.begin(), it))};
    std::size_t left{min_index - 1};
    std::size_t right{min_index};
    std::size_t const last{xs.size() - 1};
    bool do_left{min_index > 0};
    bool do_right{true};
    auto left_pos{x_axis.coord_to_pos(xs[left])};
    auto left_dist{p.x - left_pos};
    auto right_pos{x_axis.coord_to_pos(xs[right])};
    auto right_dist{right_pos - p.x};
    auto const [axis_left, axis_right] = x_axis.get_pos_range();
    double min_dist{std::numeric_limits<double>::max()};

    while (do_left || do_right)
    {
        for (; do_left && (!do_right || left_dist <= right_dist); --left)
        {
            left_pos = x_axis.coord_to_pos(xs[left]);
            left_dist = p.x - left_pos;
            auto d2{distance(left_dist, ys[left])};
            if (d2 < min_dist)
            {
                min_dist = d2;
                min_index = left;
            }
            if (left_dist > min_dist || left == 0 || left_pos < axis_left)
                do_left = false;
        }
        for (; do_right && (!do_left || right_dist <= left_dist); ++right)
        {
            right_pos = x_axis.coord_to_pos(xs[right]);
            right_dist = right_pos - p.x;
            auto d2{distance(right_dist, ys[right])};
            if (d2 < min_dist)
            {
                min_dist = d2;
                min_index = right;
            }
            if (right_dist > min_dist || right == last || right_pos > axis_right)
                do_right = false;
        }
    }
    assert(min_index < xs.size());
    return Point{xs[min_index], ys[min_index]};
}

Plotter::Plotter(Glib::RefPtr<Gtk::Application> app, Palette palette)
    : m_app(app),
      m_palette(palette),
      m_now{m_history.end()}
{
    set_can_focus(true);
    add_events(Gdk::KEY_PRESS_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK |
               Gdk::POINTER_MOTION_MASK | Gdk::SCROLL_MASK | Gdk::STRUCTURE_MASK);

    // Listen on standard input.
    m_io_connection = Glib::signal_io().connect(sigc::mem_fun(*this, &Plotter::on_read),
                                                STDIN_FILENO,
                                                Glib::IOCondition::IO_IN);
    m_io_channel = Glib::IOChannel::create_from_fd(STDIN_FILENO);
}

void Plotter::autoscale()
{
    // Return the extrema of all the points in all the vectors. Return [0.0, 1.0] if there
    // are no points.
    auto find_range = [](VV const& vv) -> std::tuple<double, double, double> {
        auto min{0.0};
        auto max{1.0};
        auto pad{0.0};
        auto found{false};
        for (auto const& v : vv)
        {
            if (v.empty())
                continue;
            auto [min_it, max_it] = std::minmax_element(v.begin(), v.end());
            min = found ? std::min(min, *min_it) : *min_it;
            max = found ? std::max(max, *max_it) : *max_it;
            pad = autoscale_padding*(std::abs(max - min));
            found = true;
        }
        // Show a span of 1.0 if there's only one point.
        if (min == max)
        {
            min -= 0.5;
            max += 0.5;
            pad = 0.0;
        }
        return {min, max, pad};
    };

    auto [x_min, x_max, x_pad] = find_range(m_xss);
    m_x_axis.set_coord_range(x_min - x_pad, x_max + x_pad);
    auto [y_min, y_max, y_pad] = find_range(m_yss);
    m_y_axis.set_coord_range(y_min - y_pad, y_max + y_pad);

    queue_draw();
}

bool Plotter::on_read(Glib::IOCondition io_cond)
{
    if ((io_cond & Glib::IOCondition::IO_IN) != Glib::IOCondition::IO_IN)
    {
        std::cerr << "Unexpected IO condition" << std::endl;
        return true;
    }

    // on_read() should only be called once.
    assert(m_xss.empty() && m_yss.empty());
    Glib::ustring data;
    auto read_state{-1};
    auto status{m_io_channel->read_to_end(data)};
    if (status != Glib::IO_STATUS_NORMAL)
        return true;

    std::istringstream is(data);
    std::string token;
    while (is)
    {
        is >> token;
        if (token == end_tag)
            break;
        if (token == data_tag)
            (++read_state % 2 == 0 ? m_xss : m_yss).push_back(std::vector<double>());
        // Ignore anything before the start tag.
        else if (read_state != -1)
            (read_state % 2 == 0 ? m_xss : m_yss).back().push_back(std::atof(token.c_str()));
    }
    m_total_points = std::accumulate(m_xss.begin(), m_xss.end(), 0U,
                                     [](std::size_t n, auto xs){ return n + xs.size(); });
    // Start with lines if there are a lot of points to avoid clutter and delay.
    m_line_style = m_total_points > big_plot_threshold ? Line_Style::lines : Line_Style::points;
    autoscale();
    record(false);

    m_io_connection.disconnect();
    close(STDIN_FILENO);
    return true;
}

void Plotter::move(double x_frac, double y_frac)
{
    if (mp_subrange)
        mp_subrange->move({x_frac*mp_subrange->width(), y_frac*mp_subrange->height()}, true);
    else
    {
        auto [x1, x2] = m_x_axis.get_pos_range();
        auto [y1, y2] = m_y_axis.get_pos_range();
        m_x_axis.move_coord_range_by_pos(x_frac*std::abs(x2 - x1));
        m_y_axis.move_coord_range_by_pos(y_frac*std::abs(y2 - y1));
    }
    queue_draw();
}

void Plotter::scale(double x_frac, double y_frac, std::optional<Point> center = std::nullopt)
{
    if (mp_subrange)
        // Scale the box instead of the range
        mp_subrange->scale(1.0/x_frac, 1.0/y_frac, center);
    else
    {
        if (center)
        {
            m_x_axis.scale_range(x_frac, center->x);
            m_y_axis.scale_range(y_frac, center->y);
        }
        else
        {
            m_x_axis.scale_range(x_frac);
            m_y_axis.scale_range(y_frac);
        }
    }
    queue_draw();
}

bool Plotter::on_key_press_event(GdkEventKey* event)
{
    auto const pan{event->state & Gdk::ModifierType::CONTROL_MASK
        ? pan_distance.second : pan_distance.first};
    auto const zoom{event->state & Gdk::ModifierType::CONTROL_MASK
        ? zoom_factor.second : zoom_factor.first};
    // Shift -> vertical zoom only
    auto const x_zoom{event->state & Gdk::ModifierType::SHIFT_MASK
        ? 1.0 : zoom};
    // Alt -> horizontal zoom only.
    auto const y_zoom{event->state & Gdk::ModifierType::MOD1_MASK
        ? 1.0 : zoom};

    switch (event->keyval)
    {
    case GDK_KEY_a:
        autoscale();
        record(false);
        break;
    case GDK_KEY_plus:
    case GDK_KEY_equal:
        scale(1.0/x_zoom, 1.0/y_zoom);
        record(true);
        break;
    case GDK_KEY_minus:
    case GDK_KEY_underscore:
        scale(x_zoom, y_zoom);
        record(true);
        break;
    case GDK_KEY_Left:
        move(-pan, 0.0);
        record(true);
        break;
    case GDK_KEY_Right:
        move(pan, 0.0);
        record(true);
        break;
    case GDK_KEY_Up:
        move(0.0, -pan);
        record(true);
        break;
    case GDK_KEY_Down:
        move(0.0, pan);
        record(true);
        break;
    case GDK_KEY_Tab:
        m_line_style = Line_Style{(static_cast<int>(m_line_style) + 1)
            % static_cast<int>(Line_Style::num)};
        queue_draw();
        break;
    case GDK_KEY_q:
    {
        auto [xlow, xhigh] = m_x_axis.get_coord_range();
        auto [ylow, yhigh] = m_y_axis.get_coord_range();
        // Send the ranges to stdout to be received by the R call.
        std::cout << m_x_axis.format(xlow, output_precision) << ' '
                  << m_x_axis.format(xhigh, output_precision) << ' '
                  << m_y_axis.format(ylow, output_precision) << ' '
                  << m_y_axis.format(yhigh, output_precision) << std::endl;
        m_app->quit();
        break;
    }
    case GDK_KEY_y:
        redo();
        break;
    case GDK_KEY_z:
        undo();
        break;
    case GDK_KEY_space:
        if (!mp_subrange)
        {
            // Enter overview mode. Note the current coordinate range before autoscaling.
            auto [x1, x2] = m_x_axis.get_coord_range();
            auto [y1, y2] = m_y_axis.get_coord_range();
            // Change the view, but don't call record(), it's only temporary.
            autoscale();
            // Set the subrange to the device coordinates of the old range on the
            // autoscaled axes.
            mp_subrange = std::make_unique<Subrange>(
                Point{m_x_axis.coord_to_pos(x1), m_y_axis.coord_to_pos(y2)},
                Point{m_x_axis.coord_to_pos(x2), m_y_axis.coord_to_pos(y1)});
        }
        else
        {
            // Leave overview mode. Set the ranges according to the subrange box.
            m_x_axis.set_coord_range_by_pos(mp_subrange->get_p1().x, mp_subrange->get_p2().x);
            m_y_axis.set_coord_range_by_pos(mp_subrange->get_p2().y, mp_subrange->get_p1().y);
            mp_subrange.reset();
            record(false);
            queue_draw();
        }
        break;
    case GDK_KEY_Escape:
        if (mp_subrange)
        {
            mp_subrange.reset();
            // Restore the range to the last recorded state.
            update(m_now);
        }
        else if (m_drag)
        {
            // Cancel any drag in progress.
            m_drag.reset();
            queue_draw();
        }
        else
        {
            // Stop displaying the closest point.
            m_closest_point.reset();
            queue_draw();
        }
        break;
    }
    return true;
}

bool Plotter::on_button_press_event(GdkEventButton* event)
{
    if (event->button != 1)
        return true;

    // If dragging starts on an axis, put the start position on the other side of the grid
    // so the full range in the other dimension is selected. Dragging on the x-axis only
    // changes the x-scale, likewise for y.
    auto [x_low, x_high] = m_x_axis.get_pos_range();
    auto [y_low, y_high] = m_y_axis.get_pos_range();
    m_drag = {{event->x < x_low ? x_high : event->x, event->y > y_low ? y_high : event->y},
              {std::max(event->x, x_low), std::min(event->y, y_low)},
              bool(event->state & Gdk::ModifierType::SHIFT_MASK)};
    // If dragging from the corner, the zoom box initially covers the whole grid. In all
    // other cases, the zoom box initially has no area.
    if (event->x < x_low && event->y > y_low)
        queue_draw();

    // Get retady to change the size or position of the subrange box.
    if (mp_subrange)
        mp_subrange->start({event->x, event->y});

    return true;
}

bool Plotter::on_motion_notify_event(GdkEventMotion* event)
{
    if (event->state & Gdk::ModifierType::BUTTON3_MASK)
        m_closest_point = find_closest_point({event->x, event->y}, m_xss[0], m_yss[0],
                                             m_x_axis, m_y_axis);
    if (m_closest_point)
        queue_draw();

    if (!m_drag || !(event->state & Gdk::ModifierType::BUTTON1_MASK))
        return true;

    auto last{m_drag->pointer};
    {
        auto [x_low, x_high] = m_x_axis.get_pos_range();
        auto [y_low, y_high] = m_y_axis.get_pos_range();
        m_drag->pointer = {clip(event->x, x_low, x_high), clip(event->y, y_high, y_low)};
        if (event->state & Gdk::ModifierType::CONTROL_MASK)
        {
            // Round off dragged distances. But if we're only changing one dimension,
            // don't round the other dimension. Rounding could cause it to move.
            if (m_x_axis.is_in_pos_range(event->x))
                m_drag->pointer.x = m_x_axis.round_pos(m_drag->pointer.x, coarse_motion_divisions);
            if (m_y_axis.is_in_pos_range(event->y))
                m_drag->pointer.y = m_y_axis.round_pos(m_drag->pointer.y, coarse_motion_divisions);
        }
    }
    auto dx{m_drag->pointer.x - last.x};
    auto dy{m_drag->pointer.y - last.y};

    if (mp_subrange)
    {
        auto last_p1{mp_subrange->get_p1()};
        auto last_p2{mp_subrange->get_p2()};
        mp_subrange->move({dx, dy});
        auto p1{mp_subrange->get_p1()};
        auto p2{mp_subrange->get_p2()};

        // Redraw the top and left borders...
        redraw(get_window(), last_p2, last_p1, p1, handle_width + 1);
        // ...then the bottom and right borders.
        redraw(get_window(), last_p1, last_p2, p2, handle_width + 1);
        // We need a little extra width to make sure the stroked rectangle is fully
        // erased.
    }
    else
        redraw(get_window(), m_drag->start, last, m_drag->pointer);

    // Push-pan can be done in both normal and overview modes.
    if (m_drag->shift)
    {
        m_x_axis.move_coord_range_by_pos(-dx);
        m_y_axis.move_coord_range_by_pos(-dy);
        queue_draw();
    }
    return true;
}

bool Plotter::on_button_release_event(GdkEventButton*)
{
    if (!mp_subrange
        && m_drag
        && m_drag->pointer.x != m_drag->start.x    // Don't resize if there wasn't motion
        && m_drag->pointer.y != m_drag->start.y)   // in both dimensions.
    {
        if (!m_drag->shift)
        {
            m_x_axis.set_coord_range_by_pos(m_drag->start.x, m_drag->pointer.x);
            m_y_axis.set_coord_range_by_pos(m_drag->start.y, m_drag->pointer.y);
        }
        record(false);
        queue_draw();
    }
    m_drag.reset();

    return true;
}

bool Plotter::on_scroll_event(GdkEventScroll* event)
{
    // Do coarse scaling with Ctrl.
    auto zoom{event->state & Gdk::ModifierType::CONTROL_MASK
        ? zoom_factor.second : zoom_factor.first};
    if (event->direction == GDK_SCROLL_UP)
        zoom = 1.0/zoom;
    // Shift+wheel scales y only.
    // Alt+wheel scales x only.
    scale(event->state & Gdk::ModifierType::SHIFT_MASK ? 1.0 : zoom,
          event->state & Gdk::ModifierType::MOD1_MASK ? 1.0 : zoom, // MOD1 is Alt.
          Point{event->x, event->y});
    record(true);
    queue_draw();
    return true;
}

bool Plotter::on_configure_event(GdkEventConfigure* event)
{
    // The window size has changed. Reset the axes and the subrange box if it's active.

    // The subrange position will change, but not the coordinates. Save the current
    // coordinates before changing the axes.
    Point p1, p2;
    if (mp_subrange)
    {
        p1 = {m_x_axis.pos_to_coord(mp_subrange->get_p1().x),
            m_y_axis.pos_to_coord(mp_subrange->get_p1().y)};
        p2 = {m_x_axis.pos_to_coord(mp_subrange->get_p2().x),
            m_y_axis.pos_to_coord(mp_subrange->get_p2().y)};
    }

    // Leave space for the range bars and the point display below the x-axis.
    m_x_axis.set_pos(grid_left_margin, event->width - border_width);
    m_y_axis.set_pos(event->height - grid_bottom_margin, border_width);

    // Set the new subrange positions from the coordinates that were saved before the axis
    // change.
    if (mp_subrange)
        mp_subrange->set({m_x_axis.coord_to_pos(p1.x), m_y_axis.coord_to_pos(p1.y)},
                         {m_x_axis.coord_to_pos(p2.x), m_y_axis.coord_to_pos(p2.y)});

    return true;
}

bool Plotter::on_draw(Context const& cr)
{
    // Clear the window
    set_color(cr, white);
    cr->rectangle(0, 0, get_width(), get_height());
    cr->fill();

    assert(m_xss.size() == m_yss.size());
    if (m_xss.empty())
        return true;

    cr->set_font_size(text_size);
    cr->select_font_face("sans", Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
    auto [x_low, x_high] = m_x_axis.get_pos_range();
    auto [y_low, y_high] = m_y_axis.get_pos_range();

    // Show the current undo frame and the total number of frames.
    {
        std::ostringstream os;
        os << std::distance(m_history.cbegin(), m_now) + 1 << '/' << m_history.size();
        set_color(cr, black);
        cr->move_to(border_width, get_height() - border_width);
        cr->show_text(os.str());
    }
    // Show the coordinates of the point closest to the pointer if available.
    if (m_closest_point
        && m_x_axis.is_in_coord_range(m_closest_point->x)
        && m_y_axis.is_in_coord_range(m_closest_point->y))
    {
        std::ostringstream os;
        os << std::setprecision(10) << *m_closest_point;
        Cairo::Matrix font_mat;
        cr->get_font_matrix(font_mat);
        cr->set_font_size(closest_point_text_size);
        // Serif font is much more readable.
        cr->select_font_face("serif", Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
        set_color(cr, black);
        Cairo::TextExtents point_ext;
        cr->get_text_extents(os.str(), point_ext);
        cr->move_to(clip(m_x_axis.coord_to_pos(m_closest_point->x) - 0.5*point_ext.width,
                         x_low, x_high - point_ext.width),
                    get_height() - border_width);
        cr->show_text(os.str());
        cr->set_font_matrix(font_mat);
        cr->select_font_face("sans", Cairo::FONT_SLANT_NORMAL, Cairo::FONT_WEIGHT_NORMAL);
    }
    // Show the axes and grid.
    {
        auto y_ticks{m_y_axis.get_ticks()};
        // Adjust the left edge to allow room for tick labels.
        auto y_label_width{0.0};
        Cairo::TextExtents ex;
        for (auto const& t : y_ticks)
        {
            if (t.label)
            {
                cr->get_text_extents(*t.label, ex);
                y_label_width = std::max(ex.width, y_label_width);
            }
        }
        // Add a small amount to make narrow gap between the labels and the range bar
        // ends. Use the difference in the font size and the text size since this is used
        // implicitly for the x-labels. See grid_bottom_margin.
        x_low = grid_left_margin + y_label_width + text_size - ex.height;
        m_x_axis.set_pos(x_low, x_high);
        // Get the x_ticks after adjusting the left endpoint to accommodate y-labels.
        auto x_ticks{m_x_axis.get_ticks()};
        // Draw the grid lines, tick marks and numbers. Do this before drawing ranges so
        // the range and its rectangle are drawn over the axis numbers.
        draw_axes_and_grid(cr, x_ticks, y_ticks, x_low, x_high, y_low, y_high);
    }

    // Draw the ranges if zoom-box dragging is in progress.
    if (mp_subrange)
        draw_ranges(cr, get_window(), mp_subrange->get_p1(), mp_subrange->get_p2(),
                    m_x_axis, m_y_axis);
    else if (m_drag)
        draw_ranges(cr, get_window(), m_drag->start, m_drag->pointer, m_x_axis, m_y_axis);

    // Mask off the area outside the graph.
    auto width{x_high - x_low};
    auto height{y_high - y_low};
    cr->rectangle(x_low, y_low, width, height);
    cr->clip();

    // Plot the data.
    std::optional<Point> closest_pos;
    if (m_closest_point)
        closest_pos = Point{m_x_axis.coord_to_pos(m_closest_point->x),
            m_y_axis.coord_to_pos(m_closest_point->y)};
    auto point_style{m_total_points > big_plot_threshold
                     ? Point_Style::square : Point_Style::round};
    for (std::size_t i{0}; i < m_xss.size(); ++i)
        draw_plot(cr, m_x_axis.coord_to_pos(m_xss[i]), m_y_axis.coord_to_pos(m_yss[i]),
                  get_color(m_palette, i, m_xss.size()), point_style, m_line_style, closest_pos);

    // Draw the interactive range box if we're in overview mode.
    if (mp_subrange)
    {
        cr->set_line_width(0.5);
        // Fill the entire rectangle...
        draw_rectangle(cr, mp_subrange->get_p1(), mp_subrange->get_p2(), zoom_box_color, true);
        // ...then draw the handles as outlines.
        for (auto side : {Side::left, Side::right, Side::top, Side::bottom})
            draw_rectangle(cr, mp_subrange->get_p1(side), mp_subrange->get_p2(side), black, false);
    }
    // Draw the zoom box if dragging is in progress.
    else if (m_drag && !m_drag->shift)
        draw_rectangle(cr, m_drag->start, m_drag->pointer, zoom_box_color, true);
    return true;
}

void Plotter::record(bool incremental)
{
    if (mp_subrange)
        return;

    // Get rid of any any redo information.
    if (m_now != m_history.end())
        m_history.erase(std::next(m_now), m_history.end());

    auto [x_min, x_max] = m_x_axis.get_coord_range();
    auto [y_min, y_max] = m_y_axis.get_coord_range();

    State new_state(x_min, x_max, y_min, y_max, incremental);
    // Don't record duplicate states. E.g. autoscale after autoscale.
    if (!m_history.empty() && new_state == m_history.back())
        return;
    // Only record the last in a string of incremental changes. E.g. keyboard or mouse
    // wheel scale changes.
    if (!m_history.empty() && m_history.back().incremental && incremental)
        m_history.back() = new_state;
    else
        m_history.push_back(new_state);
    m_now = std::prev(m_history.end());
}

void Plotter::undo()
{
    if (mp_subrange)
        return;
    assert(!m_history.empty());
    assert(m_now != m_history.end());
    if (m_now != m_history.begin())
        update(std::prev(m_now));
}

void Plotter::redo()
{
    if (mp_subrange)
        return;
    assert(!m_history.empty());
    assert(m_now != m_history.end());
    if (std::next(m_now) != m_history.end())
        update(std::next(m_now));
}

void Plotter::update(std::deque<State>::const_iterator it)
{
    m_now = it;
    m_x_axis.set_coord_range(m_now->x_min, m_now->x_max);
    m_y_axis.set_coord_range(m_now->y_min, m_now->y_max);
    queue_draw();
}

Plotter::Subrange::Subrange(Point p1, Point p2)
{
    set(p1, p2);
}

void Plotter::Subrange::set(Point p1, Point p2)
{
    m_p1 = p1;
    m_p2 = p2;
}

Point Plotter::Subrange::get_p1(Side side) const
{
    switch (side)
    {
    case Side::none:
        return m_p1;
    case Side::left:
    case Side::top:
        return {m_p1.x - handle_width, m_p1.y - handle_width};
    case Side::right:
        return {m_p2.x - handle_width, m_p1.y - handle_width};
    case Side::bottom:
        return {m_p1.x - handle_width, m_p2.y - handle_width};
    }
    return Point();
}

Point Plotter::Subrange::get_p2(Side side) const
{
    switch (side)
    {
    case Side::none:
        return m_p2;
    case Side::left:
        return {m_p1.x + handle_width, m_p2.y + handle_width};
    case Side::top:
        return {m_p2.x + handle_width, m_p1.y + handle_width};
    case Side::right:
    case Side::bottom:
        return {m_p2.x + handle_width, m_p2.y + handle_width};
    }
    return Point();
}

double Plotter::Subrange::width() const
{
    return std::abs(m_p2.x - m_p1.x);
}

double Plotter::Subrange::height() const
{
    return std::abs(m_p2.y - m_p1.y);
}

bool& Plotter::Subrange::get_side(Side side)
{
    return m_sides[static_cast<size_t>(side)];
}

void Plotter::Subrange::start(Point p)
{
    auto within = [p](Point const& p1, Point const& p2) {
        return p.x >= p1.x && p.x <= p2.x && p.y >= p1.y && p.y <= p2.y;
    };
    // Find which sides to move from the position of the pointer.
    if (within({m_p1.x + handle_width, m_p1.y + handle_width},
               {m_p2.x - handle_width, m_p2.y - handle_width}))
        m_sides.fill(true);
    else
        for (auto side : {Side::left, Side::right, Side::top, Side::bottom})
            get_side(side) = within(get_p1(side), get_p2(side));
}

void Plotter::Subrange::move(Point dp, bool all)
{
    m_p1.x += all || get_side(Side::left) ? dp.x : 0;
    m_p2.x += all || get_side(Side::right) ? dp.x : 0;
    m_p1.y += all || get_side(Side::top) ? dp.y : 0;
    m_p2.y += all || get_side(Side::bottom) ? dp.y : 0;
}

void Plotter::Subrange::scale(double x_frac, double y_frac, std::optional<Point> center)
{
    auto mid{center ? *center : Point{std::midpoint(m_p1.x, m_p2.x), std::midpoint(m_p1.y, m_p2.y)}};
    m_p1.x = x_frac*(m_p1.x - mid.x) + mid.x;
    m_p1.y = y_frac*(m_p1.y - mid.y) + mid.y;
    m_p2.x = x_frac*(m_p2.x - mid.x) + mid.x;
    m_p2.y = y_frac*(m_p2.y - mid.y) + mid.y;
}
