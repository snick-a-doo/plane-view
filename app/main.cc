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

#include <plotter.hh>

#include <gtkmm.h>

int main(int argc, char** argv)
{
    auto app = Gtk::Application::create(argc, argv, "");

    Gtk::Window window;
    Plotter plotter{app};
    plotter.plot({0, 1, 2, 3}, {1, std::sqrt(2.0), M_PI, 0.0});
    window.add(plotter);
    window.resize(500, 300);
    plotter.show();

    return app->run(window);
}
