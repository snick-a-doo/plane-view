// Copyright © 2021 Sam Varner
//
// This file is part of Plane View.
//
// 4color is free software: you can redistribute it and/or modify it under the terms of
// the GNU General Public License as published by the Free Software Foundation, either
// version 3 of the License, or (at your option) any later version.
//
// Plane View is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with 4color.
// If not, see <http://www.gnu.org/licenses/>.

#include <iomanip>
#include <iostream>
#include <cassert>
#include <numbers>
#include <numeric>
#include <sstream>

#include <plotter.hh>

auto constexpr point_radius{3.0};

struct Color
{
    int red;
    int green;
    int blue;
};

Color constexpr black{0, 0, 0};
Color constexpr light_gray{210, 210, 210};
Color constexpr white{255, 255, 255};

void set_color(Context cr, Color color)
{
    cr->set_source_rgb(color.red/255.0,
		       color.green/255.0,
		       color.blue/255.0);
}

std::pair<double, int> axis_round(double x)
{
    using namespace std;
    auto b{1 - static_cast<int>(floor(log10(x)))};
    auto lead{static_cast<int>(floor(x*pow(10, b)))};
    auto closest{100};
    auto last_l{100};
    auto prec{1};
    for (auto l : {10, 20, 25, 50, 100})
    {
	auto diff{abs(lead - l)};
	if (diff >= closest)
	{
	    prec = (last_l == 25 ? 2 : 1) - b;
	    break;
	}
	closest = diff;
	last_l = l;
    }
    std::cout << x << ' ' << b << ' ' << prec << std::endl;
    return {last_l*pow(10, -b), prec};
}

Axis::Axis()
{
}

void Axis::set_pixels(int low, int high)
{
    m_low_pos = low;
    m_high_pos = high;
}

void Axis::set_range(double low, double high)
{
    m_min = low;
    m_max = high;
}

void Axis::zoom(double factor)
{
    auto mid{std::midpoint(m_min, m_max)};
    auto dist{0.5*(m_max - m_min)/factor};
    m_min = mid - dist;
    m_max = mid + dist;
}

void Axis::set_label_pos(int pos)
{
    m_label_pos = pos;
}

Axis::VPoint Axis::ticks() const
{
    VPoint ts;
    auto [dx, x_prec] = axis_round((m_max - m_min)/m_min_ticks);
    auto low{static_cast<int>(std::ceil(m_min/dx))};
    auto high{static_cast<int>(std::floor(m_max/dx))};
    for (auto x{low}; x <= high; ++x)
    {
	std::ostringstream os;
	os << std::fixed << std::setprecision(x_prec) << x*dx;
	ts.emplace_back(to_pixels(x*dx), os.str());
    }
    return ts;
}

V Axis::to_pixels(V const& xs) const
{
    V px;
    for (auto x : xs)
	px.push_back(to_pixels(x));
    return px;
}

double Axis::to_pixels(double x) const
{
    auto scale{(m_high_pos - m_low_pos)/(m_max - m_min)};
    return m_low_pos + scale*(x - m_min);
}

void draw_grid(Context cr, Axis const& x_axis, Axis const& y_axis)
{
    auto x_ticks{x_axis.ticks()};
    auto y_ticks{y_axis.ticks()};
    for (auto x : x_ticks)
    {
	cr->move_to(x.pixel, x_axis.label_pos());
	cr->show_text(x.label);
    }
    for (auto y : y_ticks)
    {
	cr->move_to(y_axis.label_pos(), y.pixel);
	cr->show_text(y.label);
    }

    cr->rectangle(x_axis.low_pos(), y_axis.low_pos(),
		  x_axis.size(), y_axis.size());
    cr->clip();
    set_color(cr, light_gray);
    cr->rectangle(x_axis.low_pos(), y_axis.low_pos(),
		  x_axis.size(), y_axis.size());
    cr->fill();

    set_color(cr, white);
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

Plotter::Plotter()
{
    set_can_focus(true);
    add_events(Gdk::KEY_PRESS_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::STRUCTURE_MASK);
}

void Plotter::plot(V const& xs, V const& ys)
{
    m_xs = xs;
    m_ys = ys;
    auto [x_min, x_max] = std::minmax_element(m_xs.begin(), m_xs.end());
    auto [y_min, y_max] = std::minmax_element(m_ys.begin(), m_ys.end());
    m_x_axis.set_range(*x_min, *x_max);
    m_y_axis.set_range(*y_min, *y_max);
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
            Gtk::Main::quit();
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

    m_x_axis.set_pixels(label.width + 2*ex.width, event->width - ex.width);
    m_x_axis.set_label_pos(event->height - ex.width);
    m_y_axis.set_pixels(event->height - 3*ex.height, ex.height);
    m_y_axis.set_label_pos(ex.width);
    return true;
}

bool Plotter::on_draw(Context const& cr)
{
    assert(m_xs.size() == m_ys.size());
    if (m_xs.empty())
        return true;

    draw_grid(cr, m_x_axis, m_y_axis);
    draw_plot(cr, m_x_axis, m_y_axis, m_xs, m_ys, black, m_line_style);
    return true;
}
