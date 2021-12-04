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
};

Color constexpr black{0, 0, 0};
Color constexpr background{225, 225, 245};
Color constexpr white{255, 255, 255};

void set_color(Context cr, Color color)
{
    cr->set_source_rgb(color.red/255.0,
               color.green/255.0,
               color.blue/255.0);
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
        cr->move_to(px[0], py[0]);
        for (std::size_t i = 1; i < px.size(); ++i)
            cr->line_to(px[i], py[i]);
        cr->stroke();
    }
    if (style != Line_Style::lines)
    {
        cr->set_line_width(1.0);
        for (std::size_t i = 0; i < px.size(); ++i)
        {
            cr->arc(px[i], py[i], point_radius, 0.0, 2.0*std::numbers::pi);
            cr->fill();
        }
    }
}

Plotter::Plotter(Glib::RefPtr<Gtk::Application> app)
    : m_app(app)
{
    set_can_focus(true);
    add_events(Gdk::KEY_PRESS_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::STRUCTURE_MASK);

    // Can't use mkstemp() because the file is a named pipe. In general, tmpnam() is
    // unsafe because an imposter file could be created between the tmpnam() and open()
    // calls. Here, mkfifo() would fail if the file already exists
    m_pipe = std::tmpnam(nullptr);
    if (mkfifo(m_pipe.c_str(), 0600) != 0)
    {
        std::cerr << "Error creating FIFO " << m_pipe << std::endl;
        return;
    }
    auto read_fd{open(m_pipe.c_str(), O_RDWR)};
    if (read_fd == -1)
    {
        std::cerr << "Error opening " << m_pipe << std::endl;
        return;
    }
    std::cout << "Listening on " << m_pipe << std::endl;
    Glib::signal_io().connect(sigc::mem_fun(*this, &Plotter::on_read),
                              read_fd,
                              Glib::IOCondition::IO_IN);
    m_io_channel = Glib::IOChannel::create_from_fd(read_fd);
}

Plotter::~Plotter()
{
    std::cout << "unlink " << m_pipe << std::endl;
    unlink(m_pipe.c_str());
}

void Plotter::plot(V const& xs, V const& ys)
{
    m_xs = xs;
    m_ys = ys;
    autoscale();
}

void Plotter::autoscale()
{
    auto [x_min, x_max] = std::minmax_element(m_xs.begin(), m_xs.end());
    auto [y_min, y_max] = std::minmax_element(m_ys.begin(), m_ys.end());
    m_x_axis.set_range(*x_min, *x_max);
    m_y_axis.set_range(*y_min, *y_max);
}

bool Plotter::on_read(Glib::IOCondition io_cond)
{
    if ((io_cond & Glib::IOCondition::IO_IN) != Glib::IOCondition::IO_IN)
    {
        std::cerr << "Unexpected IO condition" << std::endl;
        return true;
    }

    Glib::ustring line;
    Glib::IOStatus status{Glib::IO_STATUS_NORMAL};
    while (status == Glib::IO_STATUS_NORMAL)
    {
        status = m_io_channel->read_line(line);
        std::cout << m_read_state << ' ' << line;
        if (line == "\n")
        {
            ++m_read_state;
            continue;
        }
        if (line == "start\n")
        {
            m_read_state = 0;
            continue;
        }
        if (line == "end\n")
        {
            m_xs = m_read_xs;
            m_ys = m_read_ys;
            m_read_xs.clear();
            m_read_ys.clear();
            m_read_state = -1;
            autoscale();
            queue_draw();
            return true;
        }
        switch (m_read_state)
        {
        case -1:
            break;
        case 0:
            m_read_xs.push_back(std::atof(line.c_str()));
            break;
        case 1:
            m_read_ys.push_back(std::atof(line.c_str()));
            break;
        default:
            // Multiple plots are not implemented yet.
            assert(false);
        }
    }
    return true;
}

bool Plotter::on_key_press_event(GdkEventKey* event)
{
    if (event->keyval == GDK_KEY_z)
        return true;
    else
    {
        // auto shift{event->state & Gdk::ModifierType::SHIFT_MASK};
        switch (event->keyval)
        {
        case GDK_KEY_a:
        m_zoom = 1.0;
        queue_draw();
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
        case GDK_KEY_w:
            break;
        case GDK_KEY_q:
            m_app->quit();
            break;
        default:
            return true;
        }
    }
    return true;
}

bool Plotter::on_button_press_event(GdkEventButton*)
{
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
    draw_plot(cr, m_x_axis, m_y_axis, m_xs, m_ys, black, m_line_style);
    return true;
}
