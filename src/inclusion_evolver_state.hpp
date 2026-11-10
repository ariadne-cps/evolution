#ifndef ARIADNE_EVOLUTION_INCLUSION_EVOLVER_STATE_HPP
#define ARIADNE_EVOLUTION_INCLUSION_EVOLVER_STATE_HPP
#include "evolution/differential_inclusion_evolver.hpp"
#include <algorithm>
#include <limits>
#include <stdexcept>

namespace Ariadne {
namespace EvolutionDetail {
inline Nat scheduled_step(Nat current,Nat delay) {
    if(delay>=std::numeric_limits<Nat>::digits)
        throw std::overflow_error("Inclusion schedule delay is not representable");
    Nat offset=Nat(1u)<<delay;
    if(current>std::numeric_limits<Nat>::max()-offset)
        throw std::overflow_error("Inclusion schedule step overflow");
    return current+offset;
}
}

struct ScheduledApproximator
{
    Nat step;
    InclusionIntegrator approximator;

    ScheduledApproximator(Nat s, InclusionIntegrator a) : step(s), approximator(a) {}
};

inline OutputStream& operator<<(OutputStream& os, ScheduledApproximator const& sa) {
    return os << "(" << sa.step << ":" << sa.approximator << ")"; }

struct ScheduledApproximatorComparator {
    inline bool operator() (ScheduledApproximator const& sa1, ScheduledApproximator const& sa2) {
        return (sa1.step > sa2.step); }
};

class InclusionEvolverState {
  private:
    Nat _step;
    List<ScheduledApproximator> _schedule;
    Map<InclusionIntegrator,Nat> _approximator_check_delay;
    Map<InclusionIntegrator,Nat> _approximator_global_optima_count;
    Map<InclusionIntegrator,Nat> _approximator_local_optima_count;
  public:
    InclusionEvolverState(DifferentialInclusion const& ivf, List<InputApproximation> const& approximations, IntegratorInterface const& integrator)
        : _step(0u)
    {
        InclusionIntegratorFactory factory(integrator);
        bool has_reach_rigorous_approximations = false;
        for (auto appro : approximations) {
            if (appro.is_reach_rigorous()) has_reach_rigorous_approximations = true;
            InclusionIntegrator approximator = factory.create(ivf.function(), ivf.inputs(), appro);
            _schedule.push_back(ScheduledApproximator(0u,approximator));
            _approximator_global_optima_count[approximator] = 0;
            _approximator_local_optima_count[approximator] = 0;
            _approximator_check_delay[approximator] = 0;
        }
        if (not has_reach_rigorous_approximations) {
            InclusionIntegrator approximator = factory.create(ivf.function(), ivf.inputs(), ConstantApproximation());
            _schedule.push_back(ScheduledApproximator(0u,approximator));
            _approximator_global_optima_count[approximator] = 0;
            _approximator_local_optima_count[approximator] = 0;
            _approximator_check_delay[approximator] = 0;
        }
    }

    Nat step() const { return _step; }
    List<ScheduledApproximator> const& schedule() const { return _schedule; }
    Map<InclusionIntegrator,Nat> const& check_delay() const { return _approximator_check_delay; }
    Map<InclusionIntegrator,Nat> const& global_optima_count() const { return _approximator_global_optima_count; }
    Map<InclusionIntegrator,Nat> const& local_optima_count() const { return _approximator_local_optima_count; }

    Void reset_local_optima_count() {
        for (auto entry : _approximator_local_optima_count)
            _approximator_local_optima_count[entry.first] = 0;
    }

    Void update_with_best(InclusionIntegrator const& best) {
        Nat max_reach_rigorous_index = 0;
        for (auto appro : approximators_to_use()) {
            if (appro.is_reach_rigorous()) max_reach_rigorous_index = std::max(max_reach_rigorous_index,appro.index());
        }
        for (auto appro : approximators_to_use()) {
            if (best == appro) _reset_check_delay(appro);
            else {
                if (not best.is_reach_rigorous() and appro.index() == max_reach_rigorous_index)
                    _reset_check_delay(appro);
                else _increase_check_delay(appro);
            }

            _schedule.push_back(ScheduledApproximator(
                EvolutionDetail::scheduled_step(_step,_approximator_check_delay[appro]),appro));
        }
        std::sort(_schedule.begin(),_schedule.end(),ScheduledApproximatorComparator());

        _increase_optima_count(best);
    }

    List<InclusionIntegrator> approximators_to_use() const {
        List<InclusionIntegrator> result;
        for (auto i = _schedule.size(); i > 0; --i) {
            auto entry = _schedule[i-1];
            if (entry.step == _step) {
                result.push_back(entry.approximator);
            } else {
                // No past entries survive next_step().
                break;
            }
        }
        return result;
    }

    Void next_step() {
        // Called after update_with_best: every current entry has a
        // strictly future successor, so retirement cannot empty the list.
        while (_schedule.back().step == _step) _schedule.pop_back();
        _step++;
    }

  private:
    Void _reset_check_delay(InclusionIntegrator const& approximator) { _approximator_check_delay[approximator] = 0; }
    Void _increase_check_delay(InclusionIntegrator const& approximator) { _approximator_check_delay[approximator]++; }
    Void _increase_optima_count(InclusionIntegrator const& approximator) {
        _approximator_global_optima_count[approximator]++;
        _approximator_local_optima_count[approximator]++;
    }
};

}
#endif
