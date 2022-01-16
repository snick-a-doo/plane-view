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
#include <numbers>
#include <numeric>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <plotter.hh>

/// The size of a plotted point.
auto constexpr point_radius{3.0};
/// The height of the gap below the x-axis labels.
auto constexpr margin{20};
/// The padding as fraction of the range to add to an autoscaled graph.
auto constexpr autoscale_padding{0.05};
/// The width of the adjustable border of the range box in overview mode.
auto constexpr handle_width{20};

/// An enumeration for the x and y directions.
enum class Direction {x, y};

/// RGBA color struct.
struct Color
{
    int red;
    int green;
    int blue;
    int alpha{255}; ///< Opaque by default.
};

Color constexpr black{0, 0, 0};
Color constexpr background{225, 225, 245};
Color constexpr white{255, 255, 255};
Color constexpr gray{128, 128, 128};
Color constexpr zoom_box_color{128, 128, 128, 64};

/// The colors of plotted points and lines from Brewer set 2.
std::array<Color, 8> plot_colors{
    Color{0x66, 0xc2, 0xa5}, Color{0xfc, 0x8d, 0x62},
    Color{0x8d, 0xa0, 0xcb}, Color{0xe7, 0x8a, 0xc3},
    Color{0xa6, 0xd8, 0x54}, Color{0xff, 0xd9, 0x2f},
    Color{0xe5, 0xc4, 0x94}, Color{0xb3, 0xb3, 0xb3},
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

/// Draw the axis numbers and tick marks.
static void draw_axes(Context cr, Axis const& x_axis, Axis const& y_axis)
{
    auto x_ticks{x_axis.get_ticks()};
    auto y_ticks{y_axis.get_ticks()};
    assert (!x_ticks.empty() && !y_ticks.empty());

    auto x_center = [cr](double pos, std::string& str) {
        Cairo::TextExtents ex;
        cr->get_text_extents(str, ex);
        return pos - 0.5*ex.width;
    };
    auto y_center = [cr](double pos, std::string& str) {
        Cairo::TextExtents ex;
        cr->get_text_extents(str, ex);
        return pos + 0.5*ex.height;
    };
    auto text_width = [cr](std::string& str) {
        Cairo::TextExtents ex;
        cr->get_text_extents(str, ex);
        return ex.width;
    };

    // Draw the numbers for the ticks.
    set_color(cr, black);
    for (auto x : x_ticks)
    {
        cr->move_to(x_center(x.position, x.label), x_axis.get_tick_label_pos());
        cr->show_text(x.label);
    }
    auto [x_low, x_high] = x_axis.get_pos_range();
    auto [y_low, y_high] = y_axis.get_pos_range();
    for (auto y : y_ticks)
    {
        cr->move_to(x_low - y_axis.get_tick_label_pos() - text_width(y.label),
                    y_center(y.position, y.label));
        cr->show_text(y.label);
    }

    // Draw the tick marks.
    cr->set_line_width(1.0);
    set_color(cr, black);
    for (auto x : x_ticks)
    {
        cr->move_to(x.position, y_low + 4);
        cr->line_to(x.position, y_low);
    }
    for (auto y : y_ticks)
    {
        cr->move_to(x_low - 4, y.position);
        cr->line_to(x_low, y.position);
    }
    cr->stroke();
}

/// Draw the background and grid lines.
static void draw_grid(Context cr, Axis const& x_axis, Axis const& y_axis)
{
    auto x_ticks{x_axis.get_ticks()};
    auto y_ticks{y_axis.get_ticks()};
    assert (!x_ticks.empty() && !y_ticks.empty());

    // Mask off the area outside the graph.
    auto [x_low, x_high] = x_axis.get_pos_range();
    auto [y_low, y_high] = y_axis.get_pos_range();
    auto width{x_high - x_low};
    auto height{y_high - y_low};
    cr->rectangle(x_low, y_low, width, height);
    cr->clip();

    // Draw the grid.
    set_color(cr, background);
    cr->rectangle(x_low, y_low, width, height);
    cr->fill();

    // Major grid lines
    set_color(cr, white);
    cr->set_line_width(1.0);
    for (auto x : x_ticks)
    {
        cr->move_to(x.position, y_low);
        cr->line_to(x.position, y_high);
    }
    for (auto y : y_ticks)
    {
        cr->move_to(x_low, y.position);
        cr->line_to(x_high, y.position);
    }
    cr->stroke();

    // Major grid lines
    set_color(cr, white);
    cr->set_line_width(0.5);
    if (x_ticks.size() > 1)
    {
        auto half{0.5*(x_ticks[1].position - x_ticks[0].position)};
        for (auto x : x_ticks)
        {
            cr->move_to(x.position + half, y_low);
            cr->line_to(x.position + half, y_high);
        }
    }
    if (y_ticks.size() > 1)
    {
        auto half{0.5*(y_ticks[1].position - y_ticks[0].position)};
        for (auto y : y_ticks)
        {
            cr->move_to(x_low, y.position + half);
            cr->line_to(x_high, y.position + half);
        }
    }
    cr->stroke();
}

/// Draw the data points and lines.
static void draw_plot(Context cr,
                      std::vector<double> const& pxs, std::vector<double> const& pys,
                      Color color, Line_Style style)
{
    assert(pxs.size() == pys.size());
    if (pxs.empty())
        return;

    set_color(cr, color);
    if (style != Line_Style::points)
    {
        cr->set_line_width(1.0);
        cr->move_to(pxs[0], pys[0]);
        for (std::size_t i = 1; i < pxs.size(); ++i)
            cr->line_to(pxs[i], pys[i]);
        cr->stroke();
    }
    if (style != Line_Style::lines)
        for (std::size_t i = 0; i < pxs.size(); ++i)
        {
            cr->arc(pxs[i], pys[i], point_radius, 0.0, 2.0*std::numbers::pi);
            cr->fill();
        }
}

/// Draw range bars in the margin.
void draw_range(Context cr, int x0, int y0, int x1, int y1,
                std::string const& label, Direction axis)
{
    // Draw the range lines.
    auto dx{axis == Direction::x ? 0 : margin/2};
    auto dy{axis == Direction::x ? margin/2 : 0};
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
    auto constexpr label_pad{4};
    Cairo::TextExtents label_ext;
    cr->get_text_extents(label, label_ext);
    auto label_x{axis == Direction::x
        ? std::midpoint(x0, x1) - label_ext.width/2
        : label_pad};
    auto label_y{axis == Direction::x
        ? y0 + label_ext.height/2
        : std::midpoint(y0, y1) + label_ext.height/2};
    cr->rectangle(label_x - label_pad, label_y - label_ext.height - label_pad,
                  label_ext.width + 2*label_pad, label_ext.height + 2*label_pad);
    set_color(cr, white);
    cr->fill();
    set_color(cr, black);
    cr->move_to(label_x, label_y);
    cr->show_text(label);
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

Plotter::Plotter(Glib::RefPtr<Gtk::Application> app)
    : m_app(app),
      m_now{m_history.end()}
{
    set_can_focus(true);
    add_events(Gdk::KEY_PRESS_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK |
               Gdk::BUTTON1_MOTION_MASK | Gdk::SCROLL_MASK | Gdk::STRUCTURE_MASK);

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
        auto pad{autoscale_padding};
        auto found{false};
        for (auto const& v : vv)
        {
            if (v.empty())
                continue;
            auto [min_it, max_it] = std::minmax_element(v.begin(), v.end());
            min = found ? std::min(min, *min_it) : *min_it;
            max = found ? std::max(max, *max_it) : *max_it;
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
    m_x_axis.set_coord_range(x_min, x_max, x_pad);
    auto [y_min, y_max, y_pad] = find_range(m_yss);
    m_y_axis.set_coord_range(y_min, y_max, y_pad);

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
    Glib::ustring line;
    Glib::IOStatus status{Glib::IO_STATUS_NORMAL};
    auto read_state{-1};
    while (status == Glib::IO_STATUS_NORMAL)
    {
        status = m_io_channel->read_line(line);
        if (line == "\n")
        {
            ++read_state;
            if (read_state % 2 == 0)
                m_xss.push_back(std::vector<double>());
            else
                m_yss.push_back(std::vector<double>());
            continue;
        }
        if (line == "start\n")
        {
            read_state = 0;
            m_xss.push_back(std::vector<double>());
            continue;
        }
        if (line == "end\n")
        {
            m_io_connection.disconnect();
            close(STDIN_FILENO);
            autoscale();
            record();
            return true;
        }
        if (read_state == -1)
        {
            std::cerr << "Expected 'start'\n";
            return true;
        }
        if (read_state % 2 == 0)
            m_xss.back().push_back(std::atof(line.c_str()));
        else
            m_yss.back().push_back(std::atof(line.c_str()));
    }
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
        m_x_axis.move_pos_range(x_frac*std::abs(x2 - x1));
        m_y_axis.move_pos_range(y_frac*std::abs(y2 - y1));
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
    auto shift{event->state & Gdk::ModifierType::SHIFT_MASK};

    switch (event->keyval)
    {
    case GDK_KEY_a:
        autoscale();
        record();
        break;
    case GDK_KEY_plus:
    case GDK_KEY_equal:
        scale(1.0/1.1, 1.0/1.1);
        break;
    case GDK_KEY_minus:
    case GDK_KEY_underscore:
        scale(1.1, 1.1);
        break;
    case GDK_KEY_Left:
        move(shift ? -1.0 : -0.1, 0.0);
        break;
    case GDK_KEY_Right:
        move(shift ? 1.0 : 0.1, 0.0);
        break;
    case GDK_KEY_Up:
        move(0.0, shift ? -1.0 : -0.1);
        break;
    case GDK_KEY_Down:
        move(0.0, shift ? 1.0 : 0.1);
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
        std::cout << xlow << ' ' << xhigh << ' ' << ylow << ' ' << yhigh << std::endl;
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
            m_x_axis.set_pos_range(mp_subrange->get_p1().x, mp_subrange->get_p2().x);
            m_y_axis.set_pos_range(mp_subrange->get_p2().y, mp_subrange->get_p1().y);
            mp_subrange.reset();
            record();
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
        else
        {
            // Cancel any drag in progress.
            m_drag.reset();
            queue_draw();
        }
        break;
    }
    return true;
}

bool Plotter::on_button_press_event(GdkEventButton* event)
{
    // If dragging starts on an axis, put the start position on the other side of the grid
    // so the full range in the other dimension is selected. Dragging on the x-axis only
    // changes the x-scale, likewise for y.
    auto [x_low, x_high] = m_x_axis.get_pos_range();
    auto [y_low, y_high] = m_y_axis.get_pos_range();
    m_drag = {{event->x > x_low ? event->x : x_high, event->y < y_low ? event->y : y_high},
              {event->x, event->y},
              bool(event->state & Gdk::ModifierType::SHIFT_MASK)};

    // Get retady to change the size or position of the subrange box.
    if (mp_subrange)
        mp_subrange->start({event->x, event->y});

    return true;
}

bool Plotter::on_motion_notify_event(GdkEventMotion* event)
{
    if (!m_drag)
        return true;

    auto last{m_drag->pointer};

    auto clip = [](double x, double low, double high) {
        return std::min(std::max(x, low), high);
    };
    auto [x_low, x_high] = m_x_axis.get_pos_range();
    auto [y_low, y_high] = m_y_axis.get_pos_range();
    m_drag->pointer = {clip(event->x, x_low, x_high), clip(event->y, y_high, y_low)};
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
        redraw(get_window(), last_p2, last_p1, p1, handle_width);
        // ...then the bottom and right borders.
        redraw(get_window(), last_p1, last_p2, p2, handle_width);
    }
    else
        redraw(get_window(), m_drag->start, last, m_drag->pointer);

    // Push-pan can be done in both normal and overview modes.
    if (m_drag->shift)
    {
        m_x_axis.move_pos_range(-dx);
        m_y_axis.move_pos_range(-dy);
        queue_draw();
    }
    return true;
}

bool Plotter::on_button_release_event(GdkEventButton*)
{
    // Don't resize if there wasn't any motion.
    if (!m_drag
        || m_drag->pointer.x == m_drag->start.x
        || m_drag->pointer.y == m_drag->start.y)
        return true;

    if (!mp_subrange)
    {
        if (!m_drag->shift)
        {
            m_x_axis.set_pos_range(m_drag->start.x, m_drag->pointer.x);
            m_y_axis.set_pos_range(m_drag->start.y, m_drag->pointer.y);
        }
        record();
        queue_draw();
    }
    m_drag.reset();

    return true;
}

bool Plotter::on_scroll_event(GdkEventScroll* event)
{
    auto x_frac{event->direction == GDK_SCROLL_UP ? 1.0/1.1 : 1.1};
    auto y_frac{x_frac};

    // Ctrl+wheel scales vertically only.
    if (event->state & Gdk::ModifierType::CONTROL_MASK)
        y_frac = 1.0;
    // Shift+wheel scales horizontally only.
    if (event->state & Gdk::ModifierType::SHIFT_MASK)
        x_frac = 1.0;
    // Don't record mouse wheel events for undo.
    //!! Replace last event if it's the same kind.
    scale(x_frac, y_frac, Point{event->x, event->y});
    queue_draw();
    return true;
}

bool Plotter::on_configure_event(GdkEventConfigure* event)
{
    auto cr{get_window()->create_cairo_context()};
    Cairo::TextExtents label;
    cr->get_text_extents(std::to_string(1.0), label);
    Cairo::TextExtents ex;
    cr->get_text_extents("x", ex);

    Point p1, p2;
    if (mp_subrange)
    {
        p1 = {m_x_axis.pos_to_coord(mp_subrange->get_p1().x),
            m_y_axis.pos_to_coord(mp_subrange->get_p1().y)};
        p2 = {m_x_axis.pos_to_coord(mp_subrange->get_p2().x),
            m_x_axis.pos_to_coord(mp_subrange->get_p2().y)};
    }

    m_x_axis.set_pos(label.width + 2*ex.width,
                     event->width - 1.5*ex.width,
                     event->height - 3*ex.width);
    m_y_axis.set_pos(event->height - 6*ex.height, 1.5*ex.height, ex.width);

    if (mp_subrange)
        mp_subrange->set({m_x_axis.coord_to_pos(p1.x), m_y_axis.coord_to_pos(p1.y)},
                         {m_x_axis.coord_to_pos(p2.x), m_y_axis.coord_to_pos(p2.y)});

    return true;
}

bool Plotter::on_draw(Context const& cr)
{
    set_color(cr, white);
    cr->rectangle(0, 0, get_width(), get_height());
    cr->fill();

    assert(m_xss.size() == m_yss.size());
    if (m_xss.empty())
        return true;

    // Draw the tick marks and numbers. Do this before drawing ranges so the range and its
    // rectangle are drawn over the axis numbers.
    draw_axes(cr, m_x_axis, m_y_axis);

    // Draw the ranges if zoom-box dragging is in progress.
    if (m_drag)
    {
        auto label{m_x_axis.format(std::abs(m_x_axis.pos_to_coord(m_drag->start.x)
                                            - m_x_axis.pos_to_coord(m_drag->pointer.x)), 1)};
        draw_range(cr, m_drag->start.x, get_height() - margin/2,
                   m_drag->pointer.x, get_height() - margin/2,
                   label, Direction::x);
        get_window()->invalidate_rect(
            Gdk::Rectangle(0, get_height() - margin, get_width(), margin), false);

        label = m_y_axis.format(std::abs(m_y_axis.pos_to_coord(m_drag->start.y)
                                         - m_y_axis.pos_to_coord(m_drag->pointer.y)), 1);
        draw_range(cr, margin/2, m_drag->start.y, margin/2, m_drag->pointer.y,
                   label, Direction::y);
        Cairo::TextExtents label_ext;
        cr->get_text_extents(label, label_ext);
        get_window()->invalidate_rect(
            //!! use padding constant or do inside draw_range().
            Gdk::Rectangle(0, 0, std::max<int>(label_ext.width + 6, margin), get_height()),
            false);
    }

    // Draw the grid lines. Drawing is clipped to the interior of the axes here and
    // remains in effect to the end of this function.
    draw_grid(cr, m_x_axis, m_y_axis);

    // Plot the data.
    for (std::size_t i{0}; i < m_xss.size(); ++i)
        draw_plot(cr, m_x_axis.coord_to_pos(m_xss[i]), m_y_axis.coord_to_pos(m_yss[i]),
                  plot_colors[i % plot_colors.size()], m_line_style);

    // Draw the interactive range box if we're in overview mode.
    if (mp_subrange)
    {
        // Fill the entire rectangle...
        draw_rectangle(cr, mp_subrange->get_p1(), mp_subrange->get_p2(), zoom_box_color, true);
        // ...then draw the handles with the same semi-transparent color to make the
        // handles show up darker. Note that if the same color is used for the rectangle
        // and the handles, the handles Arden't drawn when the box is the same size as the
        // grid.
        for (auto side : {Side::left, Side::right, Side::top, Side::bottom})
            draw_rectangle(cr, mp_subrange->get_p1(side), mp_subrange->get_p2(side), black, false);
    }
    // Draw the zoom box if dragging is in progress.
    else if (m_drag && !m_drag->shift)
        draw_rectangle(cr, m_drag->start, m_drag->pointer, zoom_box_color, true);
    return true;
}

void Plotter::record()
{
    // Get rid of any any redo information.
    if (m_now != m_history.end())
        m_history.erase(std::next(m_now), m_history.end());

    auto [x_min, x_max] = m_x_axis.get_coord_range();
    auto [y_min, y_max] = m_y_axis.get_coord_range();
    m_history.emplace_back(x_min, x_max, y_min, y_max);
    m_now = std::prev(m_history.end());
}

void Plotter::undo()
{
    assert(!m_history.empty());
    assert(m_now != m_history.end());
    if (m_now != m_history.begin())
        update(std::prev(m_now));
}

void Plotter::redo()
{
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
