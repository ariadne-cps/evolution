#include "evolution/hybrid_evolver.hpp"
#include "evolution/hybrid_orbit.hpp"
#include "evolution/hybrid_time.hpp"
#include "system/hybrid_automata.hpp"
#include "integration/integrator.hpp"
#include "ariadne-benchmarking.hpp"
#include <cmath>
#include <stdexcept>
#include <string>
#include <optional>
#include <vector>

using namespace Ariadne;
namespace Bench=Ariadne::Benchmarking;
namespace {
struct Options {
    double duration=80.0,step=1.0,radius=0.5,accuracy=1e-5;
    unsigned events=8u,samples=5u;
};
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
                const unsigned number=Bench::parse_unsigned<unsigned>(value);
                if(!allow_zero && number==0u) throw std::invalid_argument("expected a positive count");
                return [&setting,number] { setting=number; };
            }));
    };
    count("e","events","Maximum discrete events (default: 8)",options.events,true);
    count("n","samples","Measured runs (default: 5)",options.samples,false);
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
// Keep timing boundaries identical: only evolver.orbit is measured.
// The Runner executes warm-ups, but the separate phase clock excludes
// verification and Logging just as the original experiment did.
struct Measurement {
    double nanoseconds;
    SizeType reach,intermediate,final;
};

Measurement measure(GeneralHybridEvolver const& evolver,
                    HybridSet const& initial,HybridTime const& horizon,
                    bool describe) {
    auto const start=Bench::Clock::now();
    auto orbit=evolver.orbit(initial,horizon,Semantics::UPPER);
    auto const stop=Bench::Clock::now();
    if(orbit.final().size()==0u)
        throw std::runtime_error("evolution produced no final enclosures");
    if(describe) {
        for(auto const& enclosure:orbit.final()) {
            LOGGING_PRINTLN_AT(1,"final location="<<enclosure.location()
                <<" time="<<enclosure.time_range()
                <<" events="<<enclosure.previous_events()
                <<" bounds="<<enclosure.state_bounding_box())
        }
    }
    return {Bench::elapsed_ns(start,stop),orbit.reach().size(),
            orbit.intermediate().size(),orbit.final().size()};
}

}
int main(int argc,const char* argv[]) {
    try {
        Options options;
        Bench::CommandLine cli;
        cli.add(CommandLineContribution("evolution-watertank-benchmark",
                                        argument_parsers(options)));
        if(!cli.acquire(argc,argv))return 0;
        HybridAutomaton const tank=make_watertank();
        RealVariable height("height"),aperture("aperture");
        StringVariable valve("valve");
        HybridSet initial(DiscreteLocation(valve|"opening"),
                          {0<=height<=Real(1/16_q),0<=aperture<=Real(1/64_q)});
        HybridTime horizon(ExactDouble(options.duration),options.events);
        GeneralHybridEvolver evolver(tank);
        evolver.set_integrator(
            TaylorPicardIntegrator(StepMaximumError(options.accuracy)));
        evolver.configuration().set_maximum_step_size(
            ApproximateDouble(options.step));
        evolver.configuration().set_maximum_enclosure_radius(
            ApproximateDouble(options.radius));

        LOGGING_PRINTLN("Watertank benchmark | duration="<<options.duration
            <<" | events="<<options.events
            <<" | step="<<options.step<<" | radius="<<options.radius
            <<" | accuracy="<<options.accuracy
            <<" | concurrency="<<ThreadManager::instance().concurrency()
            <<" | samples="<<options.samples<<" | warmups="<<cli.warmups())
        LOGGING_PRINTLN("Initial: "<<initial)

        Bench::Runner runner(cli);
        Bench::Samples durations;
        std::optional<Measurement> last;
        unsigned round=0u;
        // A Runner round is one complete orbit. The sample's elapsed time is
        // recorded by the inner phase clock, not the Runner's encompassing
        // timer (which includes validation and diagnostic logging).
        auto const completed=runner.run("watertank",options.samples,
            []{},
            [&]{
                bool const measured=round>=cli.warmups();
                last=measure(evolver,initial,horizon,measured);
            },
            [&]{
                if(not last)throw std::runtime_error("missing watertank result");
                if(round>=cli.warmups()) {
                    durations.add(last->nanoseconds);
                    LOGGING_PRINTLN_AT(1,"sample="<<round-cli.warmups()+1u
                        <<" | ms="<<Bench::Table::number(last->nanoseconds/1e6)
                        <<" | reach="<<last->reach
                        <<" | intermediate="<<last->intermediate
                        <<" | final="<<last->final)
                }
                ++round;
            });
        if(completed.size()!=options.samples ||
           durations.size()!=options.samples)
            throw std::runtime_error("incorrect watertank sample count");

        Bench::Report report("Watertank orbit duration (nanoseconds)",
                             cli.summary_fields());
        report.set_parameter_header(
            "One sample = one complete orbit; timing excludes correctness "
            "checks and enclosure diagnostic output.");
        report.set_dimension_columns({
            {"Scenario"},{"Duration",0u,Bench::TableAlignment::Right},
            {"Events",0u,Bench::TableAlignment::Right},
            {"Samples",0u,Bench::TableAlignment::Right}
        });
        report.set_table_options(Bench::TableOptions::ruled());
        report.add_row({"watertank",Bench::Table::number(options.duration),
                        std::to_string(options.events),
                        std::to_string(options.samples)},
                       Bench::summarize(durations));
        report.print();
        return 0;
    } catch(std::exception const& error) {
        LOGGING_PRINTLN("benchmark_watertank: "<<error.what())
        return 1;
    }
}
