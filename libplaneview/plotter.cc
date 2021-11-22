// Copyright © 2021 Sam Varner
//
// This file is part of Plane View.
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

#include <cassert>
#include <numbers>

#include <plotter.hh>

Plotter::Plotter()
{
    set_can_focus(true);
    add_events(Gdk::KEY_PRESS_MASK | Gdk::BUTTON_PRESS_MASK | Gdk::STRUCTURE_MASK);
}

void Plotter::plot(V const& xs, V const& ys)
{
    m_xs = xs;
    m_ys = ys;
}

bool Plotter::on_key_press_event(GdkEventKey* event)
{
    if (event->keyval == GDK_KEY_z)
        return true;
    else
    {
        auto shift{event->state & Gdk::ModifierType::SHIFT_MASK};
        switch (event->keyval)
        {
        case GDK_KEY_Left:
            break;
        case GDK_KEY_Right:
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

bool Plotter::on_button_press_event(GdkEventButton* event)
{
    return true;
}

bool Plotter::on_configure_event(GdkEventConfigure* event)
{
    return true;
}

bool Plotter::on_draw(Context const& cr)
{
    assert(m_xs.size() == m_ys.size());
    if (m_xs.empty())
        return true;

    V xs(m_xs.size());
    V ys(m_ys.size());

    auto [ x_min, x_max ] = std::minmax_element(m_xs.begin(), m_xs.end());
    auto [ y_min, y_max ] = std::minmax_element(m_ys.begin(), m_ys.end());

    auto constexpr pad{0.04};
    auto x_span{*x_max - *x_min};
    auto y_span{*y_max - *y_min};
    auto x_range{(1.0 + pad)*x_span};
    auto y_range{(1.0 + pad)*y_span};
    auto x_offset{*x_min - 0.5*pad*x_span};
    auto y_offset{*y_min - 0.5*pad*y_span};
    auto x_xform = [&](double x){ return get_width() * (x - x_offset) / x_range; };
    auto y_xform = [&](double y){ return get_height() * (1.0 - (y - y_offset) / y_range); };

    std::transform(m_xs.begin(), m_xs.end(), xs.begin(), x_xform);
    std::transform(m_ys.begin(), m_ys.end(), ys.begin(), y_xform);

    auto points{m_line_style != Line_Style::lines};
    auto lines{m_line_style != Line_Style::points};

    if (lines)
    {
        cr->move_to(xs[0], ys[0]);
        for (std::size_t i = 1; i < m_xs.size(); ++i)
            cr->line_to(xs[i], ys[i]);
        cr->stroke();
    }
    if (points)
    {
        cr->set_line_width(1.0);
        for (std::size_t i = 0; i < m_xs.size(); ++i)
        {
            cr->arc(xs[i], ys[i], 2.0, 0.0, 2.0*std::numbers::pi);
            cr->stroke();
        }
    }
    return true;
}
