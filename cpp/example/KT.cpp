#include "phase.h"
#include "omnes.h"
#include "path.h"
#include "array.h"
#include "cauchy.h"
#include "gsl_interface.h"
#include "type_aliases.h"
#include "constants.h"
#include "splined_omnes.h"
#include "path_eta_pi_pi.h"
#include "enums.h"
#include "angular_average.h"
#include "dispersion_integral.h"
#include "basis_function.h"
#include "iterative_solution.h"

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>


#include<vector>
#include<string>
#include<fstream>

#include<complex>
#include<cmath>
#include<iostream>
#include<iomanip>
#include<chrono>
#include<filesystem>

    
using namespace path_eta_pi_pi;
using namespace enums;
using namespace std::complex_literals;
using type_aliases::Complex;
using namespace constants;
namespace pt = boost::property_tree;

int main(){
	pt::ptree root;
	pt::read_json("../io.json", root);
    
    // choose which process you want to calculate
    // Process decay=Process::eta_3_pi;       // eta  -> 3 pi
    // Process decay=Process::etap_3_pi;      // eta' -> 3 pi
    // Process decay=Process::etap_eta_pi_pi;   // eta' -> eta pi pi
    Process decay = Process::V_3_pi; // V -> 3 pi
    // Process decay = Process::X_3_pi; // X -> 3 pi (1-+)
    
    double cutoff=1000;
    
    double mpi{mass_pi()};
    double meta{mass_eta()};
    double metap{mass_etap()};
    double mphi{mass_phi()};
    double momega{mass_omega()};
    double mass1, mass3, mass4;
    double epsilon, epsilon_tilde, validity, validity_tilde;
    int iterations;
    int max_subs;
    double bins;

    std::filesystem::create_directory("../results/");
    
    switch (decay) {
        case Process::eta_3_pi: {
                // define input strings
            const auto phase_0=root.get<std::string>("input.eta.phase0");
            const auto phase_1=root.get<std::string>("input.eta.phase1");
            const auto phase_2=root.get<std::string>("input.eta.phase2");
            
            //define output strings
            // for eta(') --> 3pi
            const auto output_a0=root.get<std::string>("output.eta.a0");
            const auto output_b0=root.get<std::string>("output.eta.b0");
            const auto output_a1=root.get<std::string>("output.eta.a1");
            const auto output_g1=root.get<std::string>("output.eta.g1");
            const auto output_h1=root.get<std::string>("output.eta.h1");

            mass1=1.;
            mass4=meta/mpi;
            
            iterations=0;
            max_subs=0;
            
            // parameters that can be tuned to obtain a smoother curve close to thresholds
            // (choose some positive values smaller than 2):
            // matching expansion close to pseudo-threshold
            epsilon=.5;
            validity=.2;
            // matching expansion close to scattering thresholds
            epsilon_tilde=.5;
            validity_tilde=.2;

            iteration::IterationEta3Pi solution(mass1, mass4,
                                  phase_0, phase_1, phase_2,
                                  epsilon, validity, epsilon_tilde, validity_tilde,
                                  cutoff, iterations, max_subs,
                                  output_a0, output_b0, output_a1, output_g1, output_h1);

            std::ofstream myfile;
            myfile.open("../results/Omnes.dat");
            for (std::size_t i=0; i< solution.mandelstam_array.interval_s.size(); ++i) {
                myfile << solution.mandelstam_array.interval_s[i] << "\t" << std::real(solution.evaluate_omnes(solution.mandelstam_array.interval_s[i],1,Setting::above)) <<
                "\t" << std::imag(solution.evaluate_omnes(solution.mandelstam_array.interval_s[i],1,Setting::above)) << "\n";
            }
            myfile.close();

            break;
            }
        case Process::etap_3_pi: {
                // define input strings
            const auto phase_0=root.get<std::string>("input.eta.phase0");
            const auto phase_1=root.get<std::string>("input.eta.phase1");
            const auto phase_2=root.get<std::string>("input.eta.phase2");
            
            //define output strings
            // for eta(') --> 3pi
            const auto output_a0=root.get<std::string>("output.eta.a0");
            const auto output_b0=root.get<std::string>("output.eta.b0");
            const auto output_a1=root.get<std::string>("output.eta.a1");
            const auto output_g1=root.get<std::string>("output.eta.g1");
            const auto output_h1=root.get<std::string>("output.eta.h1");
            mass1=1.;
            mass4=metap/mpi;
            
            iterations=10;
            max_subs=0;
            
            // parameters that can be tuned to obtain a smoother curve close to thresholds
            // (choose some positive values smaller than 2):
            // matching expansion close to pseudo-threshold
            epsilon=.55;
            validity=.55;
            // matching expansion close to scattering thresholds
            epsilon_tilde=.2;
            validity_tilde=.4;

            iteration::IterationEtap3Pi solution(mass1, mass4,
                                  phase_0, phase_1, phase_2,
                                  epsilon, validity, epsilon_tilde, validity_tilde,
                                  cutoff, iterations, max_subs,
                                  output_a0, output_b0, output_a1, output_g1, output_h1);
            break;
            }
        case Process::etap_eta_pi_pi: {
            // define input strings
            const auto phase_0=root.get<std::string>("input.eta.phase0");
            const auto phase_1=root.get<std::string>("input.eta.phase1");
            const auto phase_2=root.get<std::string>("input.eta.phase2");
            const auto phase_eta_pi=root.get<std::string>("input.eta.phaseetapi");
            
            // for eta --> eta pi pi
            const auto output_sub_a0=root.get<std::string>("output.eta.a0_sub");
            const auto output_sub_b0=root.get<std::string>("output.eta.b0_sub");
            const auto output_sub_c0=root.get<std::string>("output.eta.c0_sub");
            const auto output_sub_d0=root.get<std::string>("output.eta.d0_sub");
            const auto output_sub_a1=root.get<std::string>("output.eta.a1_sub");
            const auto output_sub_b1=root.get<std::string>("output.eta.b1_sub");
            mass1=1.;
            mass3=meta/mpi;
            mass4=metap/mpi;
            
            iterations=3;
            max_subs=0;
            
            // parameters that can be tuned to obtain a smoother curve close to thresholds
            // (choose some positive values smaller than 2):
            // matching expansion close to pseudo-threshold
            epsilon=.2;
            validity=.4;
            /// matching expansion close to scattering thresholds
            epsilon_tilde=.2;
            validity_tilde=.4;

            iteration::IterationEtapEtaPiPi solution(mass1, mass3, mass4,
                                  phase_0, phase_1, phase_2, phase_eta_pi,
                                  epsilon, validity, epsilon_tilde, validity_tilde,
                                  cutoff, iterations, max_subs,
                                  output_sub_a0, output_sub_b0, output_sub_c0, output_sub_d0,
                                  output_sub_a1, output_sub_b1);
            break;
            }
        case Process::V_3_pi: {
                // define input strings
            const auto phase=root.get<std::string>("input.V3Pi.phase1");
            
            //define output strings
            // for V --> 3pi
            const auto output_a=root.get<std::string>("output.V3Pi.a");
            const auto output_b=root.get<std::string>("output.V3Pi.b");
            const auto output_dalitz=root.get<std::string>("output.V3Pi.dalitz");

            mass1=1.;
            // mass4=10.;

            std::vector<double> decay_masses = {3.3,3.4,3.5,3.6,3.7,3.8,3.9,4.0,4.1,4.2,4.3,4.4,4.5,4.6,4.7,4.8,4.9};
            // std::vector<double> decay_masses = {5};
            // std::generate(decay_masses.begin(),decay_masses.end(),[n = 3] () mutable {return n=0.1;});

            for (std::size_t i=0; i< decay_masses.size(); ++i) {
                std::cout << decay_masses[i] << std::endl;
            }


            for (std::size_t i=0; i< decay_masses.size(); ++i) {

                mass4=decay_masses[i];

                std::cout << decay_masses[i] << std::endl;

                std::string output_a_new;
                output_a_new = "../results/V3pi/a_" + std::to_string((double)decay_masses[i]) + ".dat";
                
                iterations=10;
                // iterations = 1;
                max_subs=1;
                
                // parameters that can be tuned to obtain a smoother curve close to thresholds
                // (choose some positive values smaller than 2):
                // matching expansion close to pseudo-threshold
                epsilon=.1;
                validity=.05;
                /// matching expansion close to scattering thresholds
                epsilon_tilde=.3;
                validity_tilde=.2;

                // bins = 5;
                //bin width for Dalitz Plot

                iteration::IterationV3Pi solution(mass1, mass4,
                                      phase,
                                      epsilon, validity, epsilon_tilde, validity_tilde,
                                      cutoff, iterations, max_subs,
                                      output_a_new, output_b);
            }

            break;
            }
        case Process::X_3_pi: {
                // define input strings
            const auto phase=root.get<std::string>("input.X3Pi.phase1");
            
            //define output strings
            // for X --> 3pi
            const auto output_a=root.get<std::string>("output.X3Pi.a");
            const auto output_b=root.get<std::string>("output.X3Pi.b");

            mass1=1.;
            mass4=std::sqrt(40.);
            // mass4=momega/mpi;
            
            // iterations=10;
            iterations = 7;
            max_subs=1;
            
            // parameters that can be tuned to obtain a smoother curve close to thresholds
            // (choose some positive values smaller than 2):
            // matching expansion close to pseudo-threshold
            epsilon=.2;
            validity=.3;
            // matching expansion close to scattering thresholds
            epsilon_tilde=.2;
            validity_tilde=.3;

            iteration::IterationX3Pi solution(mass1, mass4,
                                  phase,
                                  epsilon, validity, epsilon_tilde, validity_tilde,
                                  cutoff, iterations, max_subs,
                                  output_a, output_b);

            break;
            }
        default:
            throw("Process not implemented.");
    }
    return 0;
}
