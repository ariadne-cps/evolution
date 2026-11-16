/***************************************************************************
 *            test_enclosure_recorded_draw.cpp
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
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *  General Public License for more details.
 */

#include <cstddef>
#include <vector>

#include "../src/enclosure_recorded_draw.hpp"
#include "logging/test_runner.hpp"
#include "ariadne-testing.hpp"

using namespace Ariadne;

class OrderingCanvas final : public CanvasBase {
  public:
    std::vector<std::size_t> drawn;
    Void initialise(StringType, StringType, double, double, double, double) override {}
    Void finalise() override {}
    Void write(const char*) const override {}
    Void move_to(double x, double) override { drawn.push_back(static_cast<std::size_t>(x)); }
    Void line_to(double, double) override {}
    Void circle(double, double, double) override {}
    Void dot(double, double) override {}
    Void stroke() override {}
    Void fill() override {}
    Void set_line_width(double) override {}
    Void set_line_colour(double, double, double) override {}
    Void set_fill_colour(double, double, double) override {}
    Void set_fill_opacity(double) override {}
};

class TestEnclosureRecordedDraw : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

    void test_batch_order_and_bounded_replay() {
        // Cross two 256-element batch boundaries and include an incomplete tail.
        std::vector<std::size_t> enclosures;
        for (std::size_t i = 0; i < 600; ++i)
            enclosures.push_back(i);
        OrderingCanvas canvas;
        draw_enclosures_recorded(canvas, enclosures, [](const std::size_t& index, RecordingCanvas& recording) {
            recording.move_to(static_cast<double>(index), 0.0);
            recording.line_to(static_cast<double>(index), 1.0);
            recording.stroke();
        });
        ARIADNE_TESTING_EXPECT_EQUAL(canvas.drawn.size(), enclosures.size());
        ARIADNE_TESTING_EXPECT_EQUAL(canvas.drawn, enclosures);
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_enclosure_recorded_draw");
    ARIADNE_TESTING_ADD_TEST(TestEnclosureRecordedDraw, test_batch_order_and_bounded_replay);
    return Ariadne::Testing::run_logging_suite(ariadne_test_suite, argc, argv);
}
