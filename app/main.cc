// Copyright © 2021-2022 Sam Varner
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

#include <getopt.h>

#include <gtkmm.h>

#include <iostream>

auto constexpr info{
    "Plane View: Interactive plotter\n"
    "Version 0.2.0 © 2022 Sam Varner GPL3\n"
    "https://github.com/snick-a-doo/plane-view\n"
};

auto constexpr usage{
    "Usage: plane-view [options] [< file]\n"
    "    -h --help    Display this message and exit.\n"
    "\n"
    "plane-view is intended to be called from R. The library 'planeview' provides the\n"
    "function pview(). See its documentation for details. pview() runs a plain-view\n"
    "process and sends data to it. When plain-view is exited, pview() prints a ggplot2\n"
    "command that reproduces the plain-view plot.\n"
};

void read_options(int argc, char** argv)
{
    while (true)
    {
        static struct option options[] = {
            {"help", no_argument, nullptr, 'h'},
            {0, 0, 0, 0}};
        int index;
        auto c{getopt_long(argc, argv, "h", options, &index)};
        if (c == -1)
            break;
        switch (c)
        {
        case 0:
            return;
        case 'h':
            std::cerr << info << std::endl;
            std::cerr << usage << std::endl;
            exit(0);
        }
    }
}

int main(int argc, char** argv)
{
    read_options(argc, argv);
    auto app = Gtk::Application::create(argc, argv, "");

    Gtk::Window window;
    Plotter plotter{app};
    window.add(plotter);
    window.set_title("Plane View");
    window.resize(800, 600);
    plotter.show();

    return app->run(window);
}
