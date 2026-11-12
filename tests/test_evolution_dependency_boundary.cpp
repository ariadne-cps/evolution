/***************************************************************************
 *            test_evolution_dependency_boundary.cpp
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

/** @file tests/test_evolution_dependency_boundary.cpp
 *  @brief Test evolution dependency boundary behavior and numerical invariants.
 */

// Check that the Evolution module does not depend on the higher-level
// reasoning and dynamics modules. This is a regular Ariadne Testing suite.
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <regex>
#include <stdexcept>
#include <string>
#include <vector>

#include "ariadne-testing.hpp"

namespace fs = std::filesystem;

class TestEvolutionDependencyBoundary : public Ariadne::Testing::ContextFixture {
  public:
    using ContextFixture::ContextFixture;

  private:
    static fs::path project_root() {
        // CMake normally compiles source files using their absolute paths.
        const fs::path source_file(__FILE__);
        if (source_file.is_absolute()) {
            const auto root = source_file.parent_path().parent_path();
            if (fs::is_regular_file(root / "include/evolution/enclosure.hpp") &&
                fs::is_regular_file(root / "src/CMakeLists.txt")) {
                return root;
            }
        }
        // Support relative __FILE__ values in in-source builds, too.
        auto directory = fs::current_path();
        for (;;) {
            if (fs::is_regular_file(directory / "include/evolution/enclosure.hpp") &&
                fs::is_regular_file(directory / "src/CMakeLists.txt")) {
                return directory;
            }
            const auto parent = directory.parent_path();
            if (parent == directory || parent.empty())
                break;
            directory = parent;
        }
        throw std::runtime_error("Unable to locate the Evolution source directory");
    }

    static std::string read_file(const fs::path& file) {
        std::ifstream stream(file);
        if (!stream)
            throw std::runtime_error("Cannot open " + file.string());
        return std::string(std::istreambuf_iterator<char>(stream),
                           std::istreambuf_iterator<char>());
    }

  public:
    void test_no_upstream_reasoning_includes() {
        const auto dir = project_root() / "include/evolution";
        const std::regex forbidden(
            R"(^[ \t]*#[ \t]*include[ \t]*["<](reasoning/|dynamics/))");
        std::vector<std::string> violations;
        std::size_t scanned = 0;
        for (const auto& entry: fs::directory_iterator(dir)) {
            if (!entry.is_regular_file())
                continue;
            const auto ext = entry.path().extension();
            if (ext != ".hpp" && ext != ".cpp")
                continue;
            ++scanned;
            std::ifstream stream(entry.path());
            if (!stream)
                throw std::runtime_error("Cannot open " + entry.path().string());
            std::string line;
            std::size_t line_number = 0;
            while (std::getline(stream, line)) {
                ++line_number;
                if (std::regex_search(line, forbidden)) {
                    violations.push_back(entry.path().filename().string() + ":" +
                                         std::to_string(line_number) + ": " + line);
                }
            }
            if (stream.bad())
                throw std::runtime_error("Error reading " + entry.path().string());
        }
        for (const auto& violation: violations) {
            ARIADNE_TESTING_NOTE(violation);
        }
        ARIADNE_TESTING_PROVE(scanned > 0);
        ARIADNE_TESTING_PROVE(violations.empty());
    }

    void test_enclosure_api_does_not_expose_storage() {
        const auto source = read_file(project_root() / "include/evolution/enclosure.hpp");
        ARIADNE_TESTING_PROVE(source.find("Storage") == std::string::npos);
    }

    void test_no_reasoning_target_dependency() {
        const auto source = read_file(project_root() / "src/CMakeLists.txt");
        ARIADNE_TESTING_PROVE(source.find("ariadne-reasoning") == std::string::npos);
    }
};

int main(int argc, char** argv) {
    ARIADNE_TESTING_SUITE("test_evolution_dependency_boundary");
    ARIADNE_TESTING_ADD_TEST(TestEvolutionDependencyBoundary, test_no_upstream_reasoning_includes);
    ARIADNE_TESTING_ADD_TEST(TestEvolutionDependencyBoundary, test_enclosure_api_does_not_expose_storage);
    ARIADNE_TESTING_ADD_TEST(TestEvolutionDependencyBoundary, test_no_reasoning_target_dependency);
    return ARIADNE_TESTING_RUN_TESTS_WITH_ARGS(argc, argv);
}
