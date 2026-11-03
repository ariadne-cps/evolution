#pragma once
#include "rendering/rendering_interface.hpp"
#include "rendering/rendering_manager.hpp"
using namespace Ariadne;

class RecordingCanvas final : public CanvasInterface {
public:
    mutable SizeType writes=0u;
    mutable StringType written_filename;
    SizeType initialisations=0u;
    SizeType finalisations=0u;
    SizeType boundaries=0u;
    SizeType primitives=0u;
    SizeType styles=0u;
    StringType horizontal_axis;
    StringType vertical_axis;
    double x_lower=0.0;
    double x_upper=0.0;

    Void initialise(StringType x, StringType y,
                    double lx, double ux, double ly, double uy) override {
        static_cast<void>(ly);
        static_cast<void>(uy);
        ++initialisations;
        horizontal_axis=x;
        vertical_axis=y;
        x_lower=lx;
        x_upper=ux;
    }
    Void finalise() override { ++finalisations; }
    Void write(const char* filename) const override { ++writes; written_filename=filename; }

    Void move_to(double x, double y) override {
        static_cast<void>(x); static_cast<void>(y); ++primitives;
    }
    Void line_to(double x, double y) override {
        static_cast<void>(x); static_cast<void>(y); ++primitives;
    }
    Void circle(double x, double y, double radius) override {
        static_cast<void>(x); static_cast<void>(y);
        static_cast<void>(radius); ++primitives;
    }
    Void dot(double x, double y) override {
        static_cast<void>(x); static_cast<void>(y); ++primitives;
    }
    Void stroke() override { ++primitives; }
    Void fill() override { ++primitives; }
    Void fill_boundary(List<Point2d> const& boundary) override {
        static_cast<void>(boundary);
        ++boundaries;
    }
    Void set_line_width(double width) override {
        static_cast<void>(width); ++styles;
    }
    Void set_line_colour(double red_component, double green_component, double blue_component) override {
        static_cast<void>(red_component); static_cast<void>(green_component);
        static_cast<void>(blue_component); ++styles;
    }
    Void set_fill_opacity(double opacity) override {
        static_cast<void>(opacity); ++styles;
    }
    Void set_fill_colour(double red_component, double green_component, double blue_component) override {
        static_cast<void>(red_component); static_cast<void>(green_component);
        static_cast<void>(blue_component); ++styles;
    }
};

class RecordingBackend final : public RenderingBackendInterface {
public:
    SharedPointer<RecordingCanvas> canvas;
    mutable Nat width=0u,height=0u;
    mutable Bool animated=false;
    mutable StringType filename;
    explicit RecordingBackend(SharedPointer<RecordingCanvas> recording) : canvas(recording) { }
    SharedPointer<CanvasInterface> make_canvas(const char* name,Nat w,Nat h,Bool animation) const override {
        filename=name; width=w; height=h; animated=animation; return canvas;
    }
};

// A backend change is process-wide; restore it even if painting throws.
class BackendScope {
    RenderingBackend _previous;
public:
    explicit BackendScope(RenderingBackend backend)
        : _previous(RenderingManager::instance().backend()) {
        RenderingManager::instance().set_backend(backend);
    }
    ~BackendScope() { RenderingManager::instance().set_backend(_previous); }
};

