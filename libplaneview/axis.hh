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

#ifndef PLANE_VIEW_LIBPLANEVIEW_AXIS_HH_INCLUDED
#define PLANE_VIEW_LIBPLANEVIEW_AXIS_HH_INCLUDED

#include <optional>
#include <string>
#include <vector>

/// An axis for a graph. Handles conversion between device coordinates and graph
/// coordinates. Objects don't know if they are horizontal, vertical, or some other
/// direction.
class Axis
{
public:
    Axis() = default;
    ~Axis() = default;

    // Method and variable names use "pos" for positions in device coordinates and "coord"
    // for position in plot coordinates. Device coordinates may be pixels, but both device
    // and plot positions are doubles.

    /// Set the device positions of the ends of the axis and the tick labels. Should be
    /// called when the graph is resized. Tick label position may be on a different
    /// dimension from the endpoint positions. The client must keep track.
    void set_pos(double low_pos, double high_pos);

    /// Set the coordinates of the ends of the axis. Values are padded rounded depending
    /// on the precision of the tick labels.
    /// @param pad_fraction The fraction of the specified range to be included below and
    /// above the specified range.
    void set_coord_range(double low_coord, double high_coord, double pad_fraction = 0.0);
    /// Set the device positions of the ends of the axis. The resulting coordinate Values
    /// are rounded depending on the precision of the tick labels.
    /// @param pad_fraction The fraction of the specified range to be included below and
    /// above the specified range.
    void set_coord_range_by_pos(double low_pos, double high_pos, double pad_fraction = 0.0);
    /// Move the coordinate range by an amount that corresponds to the given position
    /// change.
    void move_coord_range_by_pos(double delta_pos);
    /// Multiply the range.
    /// @param factor The range scale factor, 0 > gives a wider range -- zooms out.
    /// @param center_pos The device position of the point that doesn't change. If not
    /// given, the center of the axis is used.
    void scale_range(double factor, std::optional<double> center_pos = std::nullopt);

    /// @return True if the passed-in device position is in the range of the positions
    /// spanned by the axis.
    bool is_in_pos_range(double pos) const;
    /// @return True if the passed-in plot coordinate is in the range of the axis.
    bool is_in_coord_range(double coord) const;
    /// @return A pair with the endpoint plot coordinates of the axis.
    std::pair<double, double> get_coord_range() const;
    /// @return A pair with the endpoint device positions of the axis.
    std::pair<double, double> get_pos_range() const;

    /// Convert from device position to plot coordinates.
    double pos_to_coord(double pos) const;
    /// Convert from plot coordinates to device position.
    double coord_to_pos(double coord) const;
    /// Convert a vector of plot coordinates to a vector of device positions3.
    std::vector<double> coord_to_pos(std::vector<double> const& coords) const;

    /// Information for rendering tick marks and their labels. The client must center the
    /// label if that's desired.
    struct Tick
    {
        double position;   ///< The device position of the tick mark.
        /// For major ticks, the formatted number to be displayed with the mark. No value
        /// for minor ticks.
        std::optional<std::string> label;
    };
    /// @return Information for rendering all the tick marks on the axis.
    std::vector<Tick> get_ticks() const;

    /// Format a number with a precision relative to the precision of the tick labels.
    /// @param x The number to format.
    /// @param extra_prec How many more digits of precision to use.
    std::string format(double x, int extra_prec = 0) const;
    double round(double pos, int extra_prec = 0) const;

private:
    /// The device position of the low end of the axis.
    double m_low_pos{0};
    /// The device position of the high end of the axis.
    double m_high_pos{100};
    /// The plot coordinate of the low end of the axis.
    double m_low_coord{0.0};
    /// The plot coordinate of the low end of the axis.
    double m_high_coord{1.0};
};

#endif // PLANE_VIEW_LIBPLANEVIEW_AXIS_HH_INCLUDED
