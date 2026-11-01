#include "pybind11.hpp"
void evolution_submodule(pybind11::module&);
void evolution_module(pybind11::module& module) {
    evolution_submodule(module);
}
