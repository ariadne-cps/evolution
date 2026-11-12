/***************************************************************************
 *            recording_canvas.hpp
 *
 *  Copyright  2026  Luca Geretti
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne Evolution.
 *
 *  Ariadne Evolution is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne Evolution is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne Evolution. If not, see <https://www.gnu.org/licenses/>.
 */

/** @file tests/recording_canvas.hpp
 *  @brief Test recording canvas and renderer backends.
 */

#pragma once
#include "rendering/rendering_interface.hpp"
#include "rendering/rendering_manager.hpp"
using namespace Ariadne;

/** @brief RecordingCanvas: test recording canvas and renderer backends. */
class RecordingCanvas final : public CanvasInterface {
  public:
    mutable SizeType writes = 0u;
    mutable StringType written_filename;
    SizeType initialisations = 0u;
    SizeType finalisations = 0u;
    SizeType boundaries = 0u;
    SizeType primitives = 0u;
    SizeType styles = 0u;
    StringType horizontal_axis;
    StringType vertical_axis;
    double x_lower = 0.0;
    double x_upper = 0.0;

    Void initialise(StringType x, StringType y,
                    double lx, double ux, double ly, double uy) override {
        static_cast<void>(ly);
        static_cast<void>(uy);
        ++initialisations;
        horizontal_axis = x;
        vertical_axis = y;
        x_lower = lx;
        x_upper = ux;
    }
    Void finalise() override { ++finalisations; }
    Void write(const char* filename) const override {
        ++writes;
        written_filename = filename;
    }

    /** @brief Perform the specified evolution operation. */
    Void move_to(double x, double y) override {
        static_cast<void>(x);
        static_cast<void>(y);
        ++primitives;
    }
    /** @brief Perform the specified evolution operation. */
    Void line_to(double x, double y) override {
        static_cast<void>(x);
        static_cast<void>(y);
        ++primitives;
    }
    /** @brief Perform the specified evolution operation. */
    Void circle(double x, double y, double radius) override {
        static_cast<void>(x);
        static_cast<void>(y);
        static_cast<void>(radius);
        ++primitives;
    }
    /** @brief Perform the specified evolution operation. */
    Void dot(double x, double y) override {
        static_cast<void>(x);
        static_cast<void>(y);
        ++primitives;
    }
    Void stroke() override { ++primitives; }
    Void fill() override { ++primitives; }
    /** @brief Perform the specified evolution operation. */
    Void fill_boundary(const List<Point2d>& boundary) override {
        static_cast<void>(boundary);
        ++boundaries;
    }
    /** @brief Perform the specified evolution operation. */
    Void set_line_width(double width) override {
        static_cast<void>(width);
        ++styles;
    }
    /** @brief Perform the specified evolution operation. */
    Void set_line_colour(double red_component, double green_component, double blue_component) override {
        static_cast<void>(red_component);
        static_cast<void>(green_component);
        static_cast<void>(blue_component);
        ++styles;
    }
    /** @brief Perform the specified evolution operation. */
    Void set_fill_opacity(double opacity) override {
        static_cast<void>(opacity);
        ++styles;
    }
    /** @brief Perform the specified evolution operation. */
    Void set_fill_colour(double red_component, double green_component, double blue_component) override {
        static_cast<void>(red_component);
        static_cast<void>(green_component);
        static_cast<void>(blue_component);
        ++styles;
    }
};

/** @brief RecordingBackend: test recording canvas and renderer backends. */
class RecordingBackend final : public RenderingBackendInterface {
  public:
    SharedPointer<RecordingCanvas> canvas;
    mutable Nat width = 0u, height = 0u;
    mutable Bool animated = false;
    mutable StringType filename;
    explicit RecordingBackend(SharedPointer<RecordingCanvas> recording) : canvas(recording) {}
    /** @brief Perform the specified evolution operation. */
    SharedPointer<CanvasInterface> make_canvas(const char* name, Nat w, Nat h, Bool animation) const override {
        filename = name;
        width = w;
        height = h;
        animated = animation;
        return canvas;
    }
};

// A backend change is process-wide; restore it even if painting throws.
class BackendScope {
    RenderingBackend _previous;

  public:
    explicit BackendScope(RenderingBackend backend)
        /** @brief Perform the specified evolution operation. */
        : _previous(RenderingManager::instance().backend()) {
        RenderingManager::instance().set_backend(backend);
    }
    ~BackendScope() { RenderingManager::instance().set_backend(_previous); }
};

