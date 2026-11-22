/***************************************************************************
 *            simulation_time.hpp
 *  Copyright 2026 Luca Geretti
 *  SPDX-License-Identifier: GPL-3.0-or-later
 ****************************************************************************/
#ifndef ARIADNE_EVOLUTION_SIMULATION_TIME_HPP
#define ARIADNE_EVOLUTION_SIMULATION_TIME_HPP
#include "algebra/vector.hpp"
namespace Ariadne::EvolutionDetail {
template<class X>
Vector<X> simulation_arguments(const Vector<X>& state, bool timed, const X& time) {
    return timed ? join(state, Vector<X>({time})) : state;
}
} // namespace Ariadne::EvolutionDetail
#endif
