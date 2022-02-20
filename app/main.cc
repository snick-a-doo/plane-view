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

#include <color.hh>
#include <plotter.hh>

#include <getopt.h>

#include <gtkmm.h>

#include <iostream>
#include <map>
#include <string>

auto constexpr info{
    "Plane View: Interactive plotter\n"
    "Version 0.2.0 © 2022 Sam Varner GPL3\n"
    "https://github.com/snick-a-doo/plane-view\n"
};

auto constexpr usage{
    "Usage: plane-view [options] [< file]\n"
    "    -h --help              Display this message and exit.\n"
    "    -p --palette <palette> Render with the specified color palette\n"
    "\n"
    "plane-view may be called from R using pview() provided by the library 'planeview'\n."
    "See its documentation for details.\n"
};

std::map<std::string, Palette> palettes{
    {"ggplot2", Palette::ggplot2},
    {"set1", Palette::brewer_set_1},
    {"set2", Palette::brewer_set_2},
    {"okabe-ito", Palette::okabe_ito},
};

void show_palettes()
{
    std::cerr << "Available palettes are\n";
    for (auto const& p : palettes)
        std::cerr << "  " << p.first << std::endl;
}

Palette read_options(int argc, char** argv)
{
    auto palette{Palette::ggplot2};
    while (true)
    {
        static struct option options[] = {
            {"palette", required_argument, nullptr, 'p'},
            {"help", no_argument, nullptr, 'h'},
            {0, 0, 0, 0}};
        int index{0};
        auto c{getopt_long(argc, argv, "p:h", options, &index)};
        if (c == -1)
            break;
        switch (c)
        {
        case 0:
            break;
        case 'p':
        {
            std::string arg{optarg};
            if (palettes.count(arg) == 0)
            {
                std::cerr << "Unknown palette '" << optarg << "'. ";
                show_palettes();
                exit(-1);
            }
            palette = palettes[arg];
            break;
        }
        case 'h':
            std::cerr << info << std::endl;
            [[fallthrough]];
        default:
            std::cerr << usage << std::endl;
            show_palettes();
            exit(0);
            break;
        }
    }
    return palette;
}

int main(int argc, char** argv)
{
    auto palette{read_options(argc, argv)};
    auto app = Gtk::Application::create();

    Gtk::Window window;
    Plotter plotter{app, palette};
    window.add(plotter);
    window.set_title("Plane View");
    window.resize(800, 600);
    plotter.show();

    return app->run(window);
}
