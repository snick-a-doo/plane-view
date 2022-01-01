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

#include <iostream>
#include <cassert>
#include <fcntl.h>
#include <fstream>
#include <numbers>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>

#include <plotter.hh>

auto constexpr point_radius{3.0};

struct Color
{
    int red;
    int green;
    int blue;
    int alpha{255};
};

Color constexpr black{0, 0, 0};
Color constexpr background{225, 225, 245};
Color constexpr white{255, 255, 255};

std::array<Color, 8> colors{
    Color{0x66, 0xc2, 0xa5}, Color{0xfc, 0x8d, 0x62},
    Color{0x8d, 0xa0, 0xcb}, Color{0xe7, 0x8a, 0xc3},
    Color{0xa6, 0xd8, 0x54}, Color{0xff, 0xd9, 0x2f},
    Color{0xe5, 0xc4, 0x94}, Color{0xb3, 0xb3, 0xb3},
};

void set_color(Context cr, Color color)
{
    cr->set_source_rgba(color.red/255.0,
                        color.green/255.0,
                        color.blue/255.0,
                        color.alpha/255.0);
}

void draw_grid(Context cr, Axis const& x_axis, Axis const& y_axis)
{
    auto x_ticks{x_axis.ticks()};
    auto y_ticks{y_axis.ticks()};
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
        cr->move_to(x_center(x.pixel, x.label), x_axis.label_pos());
        cr->show_text(x.label);
    }
    for (auto y : y_ticks)
    {
        cr->move_to(x_axis.low_pos() - y_axis.label_pos() - text_width(y.label),
                    y_center(y.pixel, y.label));
        cr->show_text(y.label);
    }

    // Draw the tick marks.
    cr->set_line_width(1.0);
    set_color(cr, black);
    for (auto x : x_ticks)
    {
        cr->move_to(x.pixel, y_axis.low_pos() + 4);
        cr->line_to(x.pixel, y_axis.low_pos());
    }
    for (auto y : y_ticks)
    {
        cr->move_to(x_axis.low_pos() - 4, y.pixel);
        cr->line_to(x_axis.low_pos(), y.pixel);
    }
    cr->stroke();

    // Mask off the area outside the graph.
    cr->rectangle(x_axis.low_pos(), y_axis.low_pos(),
                  x_axis.size(), y_axis.size());
    cr->clip();

    // Draw the grid.
    set_color(cr, background);
    cr->rectangle(x_axis.low_pos(), y_axis.low_pos(),
                  x_axis.size(), y_axis.size());
    cr->fill();

    // Major grid lines
    set_color(cr, white);
    cr->set_line_width(1.0);
    for (auto x : x_ticks)
    {
        cr->move_to(x.pixel, y_axis.low_pos());
        cr->line_to(x.pixel, y_axis.high_pos());
    }
    for (auto y : y_ticks)
    {
        cr->move_to(x_axis.low_pos(), y.pixel);
        cr->line_to(x_axis.high_pos(), y.pixel);
    }
    cr->stroke();

    // Major grid lines
    set_color(cr, white);
    cr->set_line_width(0.5);
    if (x_ticks.size() > 1)
    {
        auto half{0.5*(x_ticks[1].pixel - x_ticks[0].pixel)};
        for (auto x : x_ticks)
        {
            cr->move_to(x.pixel + half, y_axis.low_pos());
            cr->line_to(x.pixel + half, y_axis.high_pos());
        }
    }
    if (y_ticks.size() > 1)
    {
        auto half{0.5*(y_ticks[1].pixel - y_ticks[0].pixel)};
        for (auto y : y_ticks)
        {
            cr->move_to(x_axis.low_pos(), y.pixel + half);
            cr->line_to(x_axis.high_pos(), y.pixel + half);
        }
    }
    cr->stroke();
}

void draw_plot(Context cr, Axis const& x_axis, Axis const& y_axis,
           V const& xs, V const& ys, Color color, Line_Style style)
{
    assert(xs.size() == ys.size());

    set_color(cr, color);
    V px{x_axis.to_pixels(xs)};
    V py{y_axis.to_pixels(ys)};
    if (style != Line_Style::points)
    {
        cr->set_line_width(1.0);
        cr->move_to(px[0], py[0]);
        for (std::size_t i = 1; i < px.size(); ++i)
            cr->line_to(px[i], py[i]);
        cr->stroke();
    }
    if (style != Line_Style::lines)
        for (std::size_t i = 0; i < px.size(); ++i)
        {
            cr->arc(px[i], py[i], point_radius, 0.0, 2.0*std::numbers::pi);
            cr->fill();
        }
}

Plotter::Plotter(Glib::RefPtr<Gtk::Application> app)
    : m_app(app),
      m_now{m_history.end()}
{
    set_can_focus(true);
    add_events(Gdk::KEY_PRESS_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::BUTTON_RELEASE_MASK |
               Gdk::BUTTON1_MOTION_MASK | Gdk::STRUCTURE_MASK);

    // Listen on standard input.
    m_io_connection = Glib::signal_io().connect(sigc::mem_fun(*this, &Plotter::on_read),
                                                STDIN_FILENO,
                                                Glib::IOCondition::IO_IN);
    m_io_channel = Glib::IOChannel::create_from_fd(STDIN_FILENO);
}

std::pair<double, double> min_max(VV const& vv)
{
    double min = std::numeric_limits<double>::max();
    double max = std::numeric_limits<double>::min();
    for (auto const& v : vv)
    {
        auto [min_it, max_it] = std::minmax_element(v.begin(), v.end());
        min = std::min(min, *min_it);
        max = std::max(max, *max_it);
    }
    return {min, max};
}

void Plotter::autoscale()
{
    auto [x_min, x_max] = min_max(m_xs);
    m_x_axis.set_range(x_min, x_max);
    auto [y_min, y_max] = min_max(m_ys);
    m_y_axis.set_range(y_min, y_max);
    record();
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
    assert(m_xs.empty() && m_ys.empty());
    Glib::ustring line;
    Glib::IOStatus status{Glib::IO_STATUS_NORMAL};
    while (status == Glib::IO_STATUS_NORMAL)
    {
        status = m_io_channel->read_line(line);
        // std::cout << m_read_state << ' ' << line;
        if (line == "\n")
        {
            ++m_read_state;
            if (m_read_state % 2 == 0)
                m_xs.push_back(V());
            else
                m_ys.push_back(V());
            continue;
        }
        if (line == "start\n")
        {
            m_read_state = 0;
            m_xs.push_back(V());
            continue;
        }
        if (line == "end\n")
        {
            m_io_connection.disconnect();
            close(STDIN_FILENO);
            autoscale();
            return true;
        }
        if (m_read_state == -1)
        {
            std::cerr << "Expected 'start'\n";
            return true;
        }
        if (m_read_state % 2 == 0)
            m_xs.back().push_back(std::atof(line.c_str()));
        else
            m_ys.back().push_back(std::atof(line.c_str()));
    }
    return true;
}

bool Plotter::on_key_press_event(GdkEventKey* event)
{
    // auto shift{event->state & Gdk::ModifierType::SHIFT_MASK};
    switch (event->keyval)
    {
    case GDK_KEY_a:
        autoscale();
        break;
    case GDK_KEY_Left:
        m_x_axis.zoom(0.9);
        m_y_axis.zoom(0.9);
        queue_draw();
        break;
    case GDK_KEY_Right:
        m_x_axis.zoom(1.1);
        m_y_axis.zoom(1.1);
        queue_draw();
        break;
    case GDK_KEY_Up:
        break;
    case GDK_KEY_Down:
        break;
    case GDK_KEY_Page_Up:
        break;
    case GDK_KEY_Page_Down:
        break;
    case GDK_KEY_space:
        break;
    case GDK_KEY_Tab:
        m_line_style = Line_Style{(static_cast<int>(m_line_style) + 1)
            % static_cast<int>(Line_Style::num)};
        queue_draw();
        break;
    case GDK_KEY_q:
        m_app->quit();
        break;
    case GDK_KEY_y:
        redo();
        break;
    case GDK_KEY_z:
        undo();
        break;
    }
    return true;
}

bool Plotter::on_button_press_event(GdkEventButton* event)
{
    // If dragging starts on an axis, put the start position on the other side of the grid
    // so the full range in the other dimension is selected. Dragging on the x-axis only
    // changes the x-scale, likewise for y.
    m_drag_start_x = event->x > m_x_axis.low_pos() ? event->x : m_x_axis.high_pos();
    m_drag_start_y = event->y < m_y_axis.low_pos() ? event->y : m_y_axis.high_pos();

    m_drag_x = event->x;
    m_drag_y = event->y;
    return true;
}

void redraw(Glib::RefPtr<Gdk::Window> window,
            double x0, double y0, double x1, double y1, double x2, double y2)
{
    auto min3 = [](int a, int b, int c) {
        return std::min(a, std::min(b, c));
    };
    auto max3 = [](int a, int b, int c) {
        return std::max(a, std::max(b, c));
    };
    auto do_redraw = [window](int x1, int y1, int x2, int y2) {
        // Add a pixel in each direction to account for truncation.
        window->invalidate_rect(Gdk::Rectangle(x1-1, y1-1, x2-x1+2, y2-y1+2), false);
    };
    // When a corner of the zoom box is moved, we need to redraw an L-shaped region
    // that has just been added or removed from the zoom box. First, mark the vertical
    // part of the L...
    do_redraw(std::min(x1, x2), min3(y0, y1, y2), std::max(x1, x2), max3(y0, y1, y2));
    // ...then the horizontal part.
    do_redraw(min3(x0, x1, x0), std::min(y1, y2), max3(x0, x1, x2), std::max(y1, y2));
}

bool Plotter::on_motion_notify_event(GdkEventMotion* event)
{
    if (!m_drag_start_x && !m_drag_start_y)
        return true;
    auto last_x{m_drag_x};
    auto last_y{m_drag_y};
    m_drag_x = event->x;
    m_drag_y = event->y;
    redraw(get_window(),
           *m_drag_start_x, *m_drag_start_y,
           last_x, last_y,
           m_drag_x, m_drag_y);
    return true;
}

bool Plotter::on_button_release_event(GdkEventButton* event)
{
    // Don't resize if there wasn't any motion.
    if (!m_drag_start_x || m_drag_x == *m_drag_start_x
        || !m_drag_start_y || m_drag_y == *m_drag_start_y)
        return true;

    m_x_axis.set_range_pixels(*m_drag_start_x, event->x);
    m_y_axis.set_range_pixels(*m_drag_start_y, event->y);
    record();
    queue_draw();

    m_drag_start_x.reset();
    m_drag_start_y.reset();

    return true;
}

bool Plotter::on_configure_event(GdkEventConfigure* event)
{
    auto cr{get_window()->create_cairo_context()};
    Cairo::TextExtents label;
    cr->get_text_extents(std::to_string(1.0), label);
    Cairo::TextExtents ex;
    cr->get_text_extents("x", ex);

    m_x_axis.set_pixels(label.width + 2*ex.width, event->width - 1.5*ex.width);
    m_x_axis.set_label_pos(event->height - ex.width);
    m_y_axis.set_pixels(event->height - 3.5*ex.height, 1.5*ex.height);
    m_y_axis.set_label_pos(ex.width);
    return true;
}

bool Plotter::on_draw(Context const& cr)
{
    set_color(cr, white);
    cr->rectangle(0, 0, get_width(), get_height());
    cr->fill();

    assert(m_xs.size() == m_ys.size());
    if (m_xs.empty())
        return true;

    draw_grid(cr, m_x_axis, m_y_axis);
    for (std::size_t i{0}; i < m_xs.size(); ++i)
        draw_plot(cr, m_x_axis, m_y_axis,
                  m_xs[i], m_ys[i], colors[i % colors.size()], m_line_style);

    if (!m_drag_start_x && !m_drag_start_y)
        return true;

    set_color(cr, {128, 128, 128, 128});
    auto x{m_drag_start_x ? *m_drag_start_x : m_x_axis.low_pos()};
    auto y{m_drag_start_y ? *m_drag_start_y : m_y_axis.low_pos()};
    auto width{m_drag_start_x
        ? m_drag_x - *m_drag_start_x
        : m_x_axis.high_pos() - m_x_axis.low_pos()};
    auto height{m_drag_start_y
        ? m_drag_y - *m_drag_start_y
        : m_y_axis.high_pos() - m_y_axis.low_pos()};
    cr->rectangle(x, y, width, height);
    cr->fill();

    return true;
}

void Plotter::record()
{
    // Get rid of any any redo information.
    if (m_now != m_history.end())
        m_history.erase(std::next(m_now), m_history.end());

    auto [x_min, x_max] = m_x_axis.get_range();
    auto [y_min, y_max] = m_y_axis.get_range();
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
    m_x_axis.set_range(m_now->x_min, m_now->x_max);
    m_y_axis.set_range(m_now->y_min, m_now->y_max);
    queue_draw();
}
