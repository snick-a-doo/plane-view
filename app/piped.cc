// Copyright © 2021 Sam Varner
//
// This file is part of Plane View.
//
// Plane View is free software: you can redistribute it and/or modify it under the terms of
// the GNU General Public License as published by the Free Software Foundation, either
// version 3 of the License, or (at your option) any later version.
//
// 4color is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
// PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with Plane View.
// If not, see <http://www.gnu.org/licenses/>.

#include <plotter.hh>

#include <gtkmm.h>

#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

Glib::RefPtr<Glib::IOChannel> channel;

std::unique_ptr<Plotter> plotter;

bool on_read(Glib::IOCondition io_cond)
{
    std::vector<double> xs;
    std::vector<double> ys;
    int set{0};

    if ((io_cond & Glib::IOCondition::IO_IN)
	!= Glib::IOCondition::IO_IN)
    {
	std::cerr << "Unexpected IO condition" << std::endl;
	return true;
    }

    std::fstream is("/tmp/pv-pipe");
    while (is)
    {
	Glib::ustring line;
	channel->read_line(line);
	std::cout << line;
	if (line == "\n")
	{
	    ++set;
	    if (set > 1)
		break;
	    continue;
	}
	if (set == 0)
	    xs.push_back(std::atof(line.c_str()));
	else
	    ys.push_back(std::atof(line.c_str()));
    }
    std::cout << "xs" << std::endl;
    for (auto x : xs)
	std::cout << x << std::endl;
    std::cout << "ys" << std::endl;
    for (auto y : ys)
	std::cout << y << std::endl;
    plotter->plot(xs, ys);
    return true;
}

int main(int argc, char** argv)
{
    auto app = Gtk::Application::create(argc, argv, "");

    Gtk::Window window;
    plotter = std::make_unique<Plotter>();
    window.add(*plotter);
    window.resize(500, 300);
    plotter->show();

    auto fifo{"/tmp/pv-pipe"};
    if (access(fifo, F_OK) == -1)
    {
	if (mkfifo(fifo, 0666) != 0)
	{
	    std::cerr << "error creating fifo" << std::endl;
	    return -1;
	}
    }
    auto read_fd = open(fifo, O_RDWR);
    if (read_fd == -1)
    {
	std::cerr << "error opening " << fifo << std::endl;
	return -1;
    }
    Glib::signal_io().connect(sigc::ptr_fun(on_read),
			      read_fd,
			      Glib::IOCondition::IO_IN);
  channel = Glib::IOChannel::create_from_fd(read_fd);
  return app->run(window);
}
