#include "evolution/hybrid_evolver.hpp"
#include "evolution/hybrid_orbit.hpp"
#include "evolution/hybrid_time.hpp"
#include "system/hybrid_automata.hpp"
#include "integration/integrator.hpp"
#include "threading/threading_command_line_interface.hpp"
#include "threading/thread_manager.hpp"
#include "utility/command_line_interface.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using namespace Ariadne;
namespace {
struct Options {
    double duration=80.0,step=1.0,radius=0.5,accuracy=1e-5;
    unsigned events=8u,samples=5u,warmups=1u;
};
unsigned parse_count(String const& value) {
    const std::string text(value);
    if(text.empty() || text.front()=='-') throw std::invalid_argument("expected a nonnegative integer");
    std::size_t used=0u;
    auto count=std::stoull(text,&used);
    if(used!=text.size() || count>std::numeric_limits<unsigned>::max()) throw std::invalid_argument("integer out of range");
    return static_cast<unsigned>(count);
}
List<ArgumentParser> argument_parsers(Options& options) {
    List<ArgumentParser> parsers;
    auto positive=[&](String const& short_id,String const& name,String const& description,double& setting) {
        parsers.append(make_valued_argument_parser(short_id,name,description,
            [&setting](String const& value)->VoidFunction {
                std::size_t used=0u;
                const std::string text(value);
                const double number=std::stod(text,&used);
                if(used!=text.size() || !std::isfinite(number) || number<=0.0) throw std::invalid_argument("expected a finite positive number");
                return [&setting,number] { setting=number; };
            }));
    };
    positive("d","duration","Continuous time horizon (default: 80)",options.duration);
    positive("p","step","Maximum integration step (default: 1)",options.step);
    positive("r","radius","Maximum enclosure radius (default: 0.5)",options.radius);
    positive("a","accuracy","Integrator step maximum error (default: 1e-5)",options.accuracy);
    auto count=[&](String const& short_id,String const& name,String const& description,unsigned& setting,bool allow_zero) {
        parsers.append(make_valued_argument_parser(short_id,name,description,
            [&setting,allow_zero](String const& value)->VoidFunction {
                const unsigned number=parse_count(value);
                if(!allow_zero && number==0u) throw std::invalid_argument("expected a positive count");
                return [&setting,number] { setting=number; };
            }));
    };
    count("e","events","Maximum discrete events (default: 8)",options.events,true);
    count("n","samples","Measured runs (default: 5)",options.samples,false);
    count("w","warmups","Unmeasured warmup runs (default: 1)",options.warmups,true);
    return parsers;
}
HybridAutomaton make_watertank() {
    RealVariable height("height"),aperture("aperture");
    StringVariable valve("valve");
    DiscreteLocation open(valve|"open"),opening(valve|"opening"),closed(valve|"closed"),closing(valve|"closing");
    HybridAutomaton tank("watertank");
    RealExpression flow=-Real(1/50_q)*height+Real(3/10_q)*aperture;
    tank.new_mode(open,{dot(height)=flow,dot(aperture)=0});
    tank.new_mode(closed,{dot(height)=flow,dot(aperture)=0});
    tank.new_mode(opening,{dot(height)=flow,dot(aperture)=Real(1/4_q)});
    tank.new_mode(closing,{dot(height)=flow,dot(aperture)=-Real(1/4_q)});
    tank.new_transition(closed,DiscreteEvent("start_opening"),opening,{next(height)=height,next(aperture)=aperture},height<=Real(11/2_q),EventKind::URGENT);
    tank.new_transition(open,DiscreteEvent("start_closing"),closing,{next(height)=height,next(aperture)=aperture},height>=8,EventKind::URGENT);
    tank.new_transition(opening,DiscreteEvent("finished_opening"),open,{next(height)=height,next(aperture)=aperture},aperture>=1,EventKind::URGENT);
    tank.new_transition(closing,DiscreteEvent("finished_closing"),closed,{next(height)=height,next(aperture)=aperture},aperture<=0,EventKind::URGENT);
    return tank;
}
struct Measurement { double milliseconds; SizeType reach,intermediate,final; };
Measurement measure(GeneralHybridEvolver const& evolver,HybridSet const& initial,HybridTime const& horizon,bool measured) {
    const auto start=std::chrono::steady_clock::now();
    auto orbit=evolver.orbit(initial,horizon,Semantics::UPPER);
    const auto stop=std::chrono::steady_clock::now();
    if(orbit.final().size()==0u) throw std::runtime_error("evolution produced no final enclosures");
    if(measured) {
        for(auto const& enclosure : orbit.final()) {
            LOGGING_PRINTLN_AT(1,"final location="<<enclosure.location()<<" time="<<enclosure.time_range()
                <<" events="<<enclosure.previous_events()<<" bounds="<<enclosure.state_bounding_box())
        }
    }
    return {std::chrono::duration<double,std::milli>(stop-start).count(),orbit.reach().size(),orbit.intermediate().size(),orbit.final().size()};
}
}
int main(int argc,const char* argv[]) {
    ThreadManager::instance();
    Logger::instance().configuration().set_verbosity(1u);
    try {
        Options options;
        ThreadingCommandLineInterface cli;
        cli.add(CommandLineContribution("evolution-watertank-benchmark",argument_parsers(options)));
        if(!cli.acquire(argc,argv)) {
            for(int i=1;i<argc;++i) { if(std::string(argv[i])=="--help" || std::string(argv[i])=="-h") return 0; }
            return 1;
        }
        auto tank=make_watertank();
        RealVariable height("height"),aperture("aperture");
        StringVariable valve("valve");
        HybridSet initial(DiscreteLocation(valve|"opening"),{0<=height<=Real(1/16_q),0<=aperture<=Real(1/64_q)});
        HybridTime horizon(ExactDouble(options.duration),options.events);
        GeneralHybridEvolver evolver(tank);
        evolver.set_integrator(TaylorPicardIntegrator(StepMaximumError(options.accuracy)));
        evolver.configuration().set_maximum_step_size(ApproximateDouble(options.step));
        evolver.configuration().set_maximum_enclosure_radius(ApproximateDouble(options.radius));
        LOGGING_PRINTLN("Watertank benchmark | duration="<<options.duration<<" | events="<<options.events
            <<" | step="<<options.step<<" | radius="<<options.radius<<" | accuracy="<<options.accuracy
            <<" | concurrency="<<ThreadManager::instance().concurrency()<<" | samples="<<options.samples<<" | warmups="<<options.warmups)
        LOGGING_PRINTLN("Initial: "<<initial)
        for(unsigned run=0u;run<options.warmups;++run) { measure(evolver,initial,horizon,false); }
        std::vector<double> durations;
        for(unsigned run=0u;run<options.samples;++run) {
            auto result=measure(evolver,initial,horizon,true);
            durations.push_back(result.milliseconds);
            LOGGING_PRINTLN("sample="<<run+1u<<" | ms="<<std::fixed<<std::setprecision(3)<<result.milliseconds
                <<" | reach="<<result.reach<<" | intermediate="<<result.intermediate<<" | final="<<result.final)
        }
        std::sort(durations.begin(),durations.end());
        const auto middle=durations.size()/2u;
        const double median=durations.size()%2u ? durations[middle] : (durations[middle-1u]+durations[middle])/2.0;
        LOGGING_PRINTLN("Watertank timing (ms) | median="<<std::fixed<<std::setprecision(3)<<median
            <<" | min="<<durations.front()<<" | max="<<durations.back())
        return 0;
    } catch(std::exception const& error) {
        LOGGING_PRINTLN("benchmark_watertank: "<<error.what())
        return 1;
    }
}
