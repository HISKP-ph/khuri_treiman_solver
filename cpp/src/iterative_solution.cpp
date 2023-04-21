#include "iterative_solution.h"

namespace iteration {

Iteration::Iteration(double mass_1, double mass_2,
          double mass_3, double mass_4,
          double epsilon, double validity,
          double epsilon_tilde, double validity_tilde,
          double cutoff, int iterations, int max_subs)
:
mass_1{mass_1}, mass_2{mass_2},
mass_3{mass_3}, decay_mass{mass_4},
epsilon{epsilon}, validity{validity},
epsilon_tilde{epsilon_tilde}, validity_tilde{validity_tilde},
cutoff{cutoff}, iterations{iterations}, max_subs{max_subs}
{}

IterationEta3Pi::IterationEta3Pi(double mass_1, double mass_4,
          const std::string& phase_0,
          const std::string& phase_1,
          const std::string& phase_2,
          double epsilon, double validity,
          double epsilon_tilde, double validity_tilde,
          double cutoff, int iterations, int max_subs,
          const std::string& output_a0,
          const std::string& output_b0,
          const std::string& output_a1,
          const std::string& output_g1,
          const std::string& output_h1)
:
Iteration(mass_1,mass_1,mass_1,mass_4,epsilon,validity,epsilon_tilde,validity_tilde,cutoff,iterations,max_subs),
output_a0{output_a0}, output_b0{output_b0},
output_a1{output_a1}, output_g1{output_g1}, output_h1{output_h1},
mandelstam_array(mass_1,mass_4,cutoff),
phases(phase_0,phase_1,phase_2),
omnes(mass_1,decay_mass, cutoff, phases)
{solve_khuri_treiman();}

IterationEtap3Pi::IterationEtap3Pi(double mass_1, double mass_4,
          const std::string& phase_0,
          const std::string& phase_1,
          const std::string& phase_2,
          double epsilon, double validity,
          double epsilon_tilde, double validity_tilde,
          double cutoff, int iterations, int max_subs,
          const std::string& output_a0,
          const std::string& output_b0,
          const std::string& output_a1,
          const std::string& output_g1,
          const std::string& output_h1)
:
Iteration(mass_1,mass_1,mass_1,mass_4,epsilon,validity,epsilon_tilde,validity_tilde,cutoff,iterations,max_subs),
output_a0{output_a0}, output_b0{output_b0},
output_a1{output_a1}, output_g1{output_g1}, output_h1{output_h1},
mandelstam_array(mass_1,mass_4,cutoff),
phases(phase_0,phase_1,phase_2),
omnes(mass_1,decay_mass, cutoff, phases)
{solve_khuri_treiman();}

IterationEtapEtaPiPi::IterationEtapEtaPiPi(double mass_1,
          double mass_3, double mass_4,
          const std::string& phase_0,
          const std::string& phase_1,
          const std::string& phase_2,
          const std::string& phase_eta_pi,
          double epsilon, double validity,
          double epsilon_tilde, double validity_tilde,
          double cutoff, int iterations, int max_subs,
          const std::string& output_sub_a0,
          const std::string& output_sub_b0,
          const std::string& output_sub_c0,
          const std::string& output_sub_d0,
          const std::string& output_sub_a1,
          const std::string& output_sub_b1)
:
Iteration(mass_1,mass_1,mass_3,mass_4,epsilon,validity,epsilon_tilde,validity_tilde,cutoff,iterations,max_subs),
output_sub_a0{output_sub_a0}, output_sub_b0{output_sub_b0}, output_sub_c0{output_sub_c0},
output_sub_d0{output_sub_d0}, output_sub_a1{output_sub_a1}, output_sub_b1{output_sub_b1},
mandelstam_array(mass_1,mass_3,mass_4,cutoff),
phases(phase_0,phase_1,phase_2,phase_eta_pi),
omnes(mass_1,mass_3,decay_mass, cutoff, phases)
{solve_khuri_treiman();}

IterationV3Pi::IterationV3Pi(double mass_1, double mass_4,
          const std::string& phase,
          double epsilon, double validity,
          double epsilon_tilde, double validity_tilde,
          double cutoff, int iterations, int max_subs,
          const std::string& output_a,
          const std::string& output_b)
:
Iteration(mass_1,mass_1,mass_1,mass_4,epsilon,validity,epsilon_tilde,validity_tilde,cutoff,iterations,max_subs),
output_a{output_a}, output_b{output_b},
mandelstam_array(mass_1,mass_4,cutoff),
phases(phase),
omnes(mass_1,decay_mass, cutoff, phases)
{solve_khuri_treiman();}

IterationX3Pi::IterationX3Pi(double mass_1, double mass_4,
          const std::string& phase,
          double epsilon, double validity,
          double epsilon_tilde, double validity_tilde,
          double cutoff, int iterations, int max_subs,
          const std::string& output_a,
          const std::string& output_b)
:
Iteration(mass_1,mass_1,mass_1,mass_4,epsilon,validity,epsilon_tilde,validity_tilde,cutoff,iterations,max_subs),
output_a{output_a}, output_b{output_b},
mandelstam_array(mass_1,mass_4,cutoff),
phases(phase),
omnes(mass_1,decay_mass, cutoff, phases)
{solve_khuri_treiman();}

IterationT3Pi::IterationT3Pi(double mass_1, double mass_4,
          const std::string& phase,
          double epsilon, double validity,
          double epsilon_tilde, double validity_tilde,
          double cutoff, int iterations, int max_subs,
          const std::string& output_a,
          const std::string& output_b)
:
Iteration(mass_1,mass_1,mass_1,mass_4,epsilon,validity,epsilon_tilde,validity_tilde,cutoff,iterations,max_subs),
output_a{output_a}, output_b{output_b},
mandelstam_array(mass_1,mass_4,cutoff),
phases(phase),
omnes(mass_1,decay_mass, cutoff, phases)
{solve_khuri_treiman();}

void IterationEta3Pi::write_output(SubtractionConstant sub, int iteration,
                        const std::string output_file){
    std::ofstream myfile;
    myfile.open(output_file);
    switch (sub) {
        case SubtractionConstant::a0: case SubtractionConstant::b0: case SubtractionConstant::a1:
            for (std::size_t i=0; i< mandelstam_array.interval_s.size(); ++i) {
            myfile << std::setprecision(15)<<mandelstam_array.interval_s[i] << "\t" <<
            std::real(basis_list[iteration](mandelstam_array.interval_s[i], 0, sub, Setting::above))<< "\t"<<
            std::imag(basis_list[iteration](mandelstam_array.interval_s[i], 0, sub, Setting::above))<< "\t"<<
            std::real(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t"<<
            std::imag(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t"<<
            std::real(basis_list[iteration](mandelstam_array.interval_s[i], 2, sub, Setting::above))<< "\t"<<
            std::imag(basis_list[iteration](mandelstam_array.interval_s[i], 2, sub, Setting::above))<< "\t"<<"\n";
            }
            break;
            
            
        case SubtractionConstant::g1:
            for (std::size_t i=0; i< mandelstam_array.interval_s.size(); ++i) {
            myfile << std::setprecision(15)<<mandelstam_array.interval_s[i] << "\t" <<
            0.<< "\t"<<
            0.<< "\t"<<
            std::real(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t"<<
            std::imag(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t"<<
            0.<< "\t"<<
            0.<< "\t"<<"\n";
            }
            break;
            
            
        case SubtractionConstant::h1:
            for (std::size_t i=0; i< mandelstam_array.interval_s.size(); ++i) {
            myfile << std::setprecision(15)<<mandelstam_array.interval_s[i] << "\t" <<
                0.<< "\t"<<
                0.<< "\t"<<
            std::real(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t"<<
            std::imag(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t"<<
            std::real(basis_list[iteration](mandelstam_array.interval_s[i], 2, sub, Setting::above))<< "\t"<<
            std::imag(basis_list[iteration](mandelstam_array.interval_s[i], 2, sub, Setting::above))<< "\t"<<"\n";
            }
            break;
            
        default:
            throw std::domain_error{"Subtraction constant does not belong to desired decay!"};
            break;
    }
    myfile.close();
    return;
}

void IterationEtap3Pi::write_output(SubtractionConstant sub, int iteration,
                        const std::string output_file){
    std::ofstream myfile;
    myfile.open(output_file);
    switch (sub) {
        case SubtractionConstant::a0: case SubtractionConstant::b0: case SubtractionConstant::a1:
            for (std::size_t i=0; i< mandelstam_array.interval_s.size(); ++i) {
            myfile << std::setprecision(15)<<mandelstam_array.interval_s[i] << "\t" <<
            std::real(basis_list[iteration](mandelstam_array.interval_s[i], 0, sub, Setting::above))<< "\t"<<
            std::imag(basis_list[iteration](mandelstam_array.interval_s[i], 0, sub, Setting::above))<< "\t"<<
            std::real(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t"<<
            std::imag(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t"<<
            std::real(basis_list[iteration](mandelstam_array.interval_s[i], 2, sub, Setting::above))<< "\t"<<
            std::imag(basis_list[iteration](mandelstam_array.interval_s[i], 2, sub, Setting::above))<< "\t"<<"\n";
            }
            break;
            
            
        case SubtractionConstant::g1:
            for (std::size_t i=0; i< mandelstam_array.interval_s.size(); ++i) {
            myfile << std::setprecision(15)<<mandelstam_array.interval_s[i] << "\t" <<
            0.<< "\t"<<
            0.<< "\t"<<
            std::real(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t"<<
            std::imag(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t"<<
            0.<< "\t"<<
            0.<< "\t"<<"\n";
            }
            break;
            
            
        case SubtractionConstant::h1:
            for (std::size_t i=0; i< mandelstam_array.interval_s.size(); ++i) {
            myfile << std::setprecision(15)<<mandelstam_array.interval_s[i] << "\t" <<
            0.<< "\t"<<
            0.<< "\t"<<
            std::real(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t"<<
            std::imag(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t"<<
            std::real(basis_list[iteration](mandelstam_array.interval_s[i], 2, sub, Setting::above))<< "\t"<<
            std::imag(basis_list[iteration](mandelstam_array.interval_s[i], 2, sub, Setting::above))<< "\t"<<"\n";
            }
            break;
            
        default:
            throw std::domain_error{"Subtraction constant does not belong to desired decay!"};
            break;
    }
    myfile.close();
    return;
}

void IterationEtapEtaPiPi::write_output(SubtractionConstant sub, int iteration,
                        const std::string output_file){
    std::ofstream myfile;
    myfile.open(output_file);
    switch (sub) {
        case SubtractionConstant::sub_a0: case SubtractionConstant::sub_b0: case SubtractionConstant::sub_c0: case SubtractionConstant::sub_d0:
            for (std::size_t i=0; i< mandelstam_array.interval_s.size(); ++i) {
            myfile << std::setprecision(15)<<mandelstam_array.interval_s[i] << "\t" <<
            real(basis_list[iteration](mandelstam_array.interval_s[i], 000, sub, Setting::above))<< "\t"<<
            imag(basis_list[iteration](mandelstam_array.interval_s[i], 000, sub, Setting::above))<< "\t"<<
            real(basis_list[iteration](mandelstam_array.interval_s[i], 010, sub, Setting::above))<< "\t"<<
            imag(basis_list[iteration](mandelstam_array.interval_s[i], 010, sub, Setting::above))<< "\t"<< "\t"<<"\n";
            }
            break;
            
        case SubtractionConstant::sub_a1: case SubtractionConstant::sub_b1:
            for (std::size_t i=0; i< mandelstam_array.interval_s.size(); ++i) {
            myfile << std::setprecision(15)<<mandelstam_array.interval_s[i] << "\t" <<
            real(basis_list[iteration](mandelstam_array.interval_s[i], 111, sub, Setting::above))<< "\t"<<
            imag(basis_list[iteration](mandelstam_array.interval_s[i], 111, sub, Setting::above))<< "\t"<<
            real(basis_list[iteration](mandelstam_array.interval_s[i], 110, sub, Setting::above))<< "\t"<<
            imag(basis_list[iteration](mandelstam_array.interval_s[i], 110, sub, Setting::above))<< "\t"<< "\t"<<"\n";
            }
            break;
            
        default:
            throw std::domain_error{"Subtraction constant does not belong to desired decay!"};
            break;
    }
    myfile.close();
    return;
}

void IterationV3Pi::write_output(SubtractionConstant sub, int iteration,
                        const std::string output_file){
    std::ofstream myfile;
    myfile.open(output_file);
    for (std::size_t i=0; i< mandelstam_array.interval_s.size(); ++i) {
        myfile << std::setprecision(15)<<mandelstam_array.interval_s[i] << "\t" <<
        std::real(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t" <<
        std::imag(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\n";
    }
    myfile.close();
    return;
}

void IterationX3Pi::write_output(SubtractionConstant sub, int iteration,
                        const std::string output_file){
    std::ofstream myfile;
    myfile.open(output_file);
    for (std::size_t i=0; i< mandelstam_array.interval_s.size(); ++i) {
        myfile << std::setprecision(15)<<mandelstam_array.interval_s[i] << "\t" <<
        std::real(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t"<<
        std::imag(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t"<< "\n";
    }
    myfile.close();
    return;
}

void IterationT3Pi::write_output(SubtractionConstant sub, int iteration,
                        const std::string output_file){
    std::ofstream myfile;
    myfile.open(output_file);
    for (std::size_t i=0; i< mandelstam_array.interval_s.size(); ++i) {
        myfile << std::setprecision(15)<<mandelstam_array.interval_s[i] << "\t" <<
        std::real(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t"<<
        std::imag(basis_list[iteration](mandelstam_array.interval_s[i], 1, sub, Setting::above))<< "\t"<< "\n";
    }
    myfile.close();
    return;
}

Complex IterationEta3Pi::evaluate_omnes(double s, int isospin, Setting evaluation){
    return omnes(s,isospin,evaluation);
}

Complex IterationEtap3Pi::evaluate_omnes(double s, int isospin, Setting evaluation){
    return omnes(s,isospin,evaluation);
}

Complex IterationEtapEtaPiPi::evaluate_omnes(double s, int isospin, Setting evaluation){
    return omnes(s,isospin,evaluation);
}

Complex IterationV3Pi::evaluate_omnes(double s, int isospin, Setting evaluation){
    return omnes(s,isospin,evaluation);
}

Complex IterationX3Pi::evaluate_omnes(double s, int isospin, Setting evaluation){
    return omnes(s,isospin,evaluation);
}

Complex IterationT3Pi::evaluate_omnes(double s, int isospin, Setting evaluation){
    return omnes(s,isospin,evaluation);
}

void IterationEta3Pi::solve_khuri_treiman(){
    if ((max_subs != 0) ){
        throw std::domain_error{"Eta -> 3Pi is only implemented for a constant subtraction scheme. Set max_subs to zero."};
    }
    auto start1 = std::chrono::high_resolution_clock::now();
    
    const auto omn
    {[this](double s, int isospin, Setting evaluation){
        return evaluate_omnes(s,isospin,evaluation);
    }};    

    auto starthom = std::chrono::high_resolution_clock::now();
    std::cout<<"Building Homogeneous solution ...\n";

    basis::HomogeneousSolutionEta3Pi homsol(mass_1,decay_mass, omn);
    
    const auto hom{[homsol](double s, int isospin, SubtractionConstant sub, Setting evaluation){
        return homsol.initial_guess(s, isospin, sub, evaluation);}};
    
    basis_list.push_back(hom);
    
    auto finishhom = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedhom = finishhom - starthom;
    std::cout<<"Finished in "<< elapsedhom.count()  << " seconds!" <<'\n';
    
    for (int j=0; j<iterations; j++) {
        std::cout<<"----------- Iteration "<< j+1 << " ------------"<<"\n";
        //-- Inhomogeneity ---------------------------------------------
        auto start2 = std::chrono::high_resolution_clock::now();
        std::cout<<"Initializing BasisTilde ...\n";
       
        // Initialize the correct basisfunction
        basis::BasisTildeEta3Pi tild(mass_1,decay_mass,
                                        basis_list[j], cutoff, j,
                                        epsilon_tilde, validity_tilde);
                
        const auto matched_tilde{[tild](double s, int isospin, SubtractionConstant sub){
            return tild.matched_tilde(s, isospin, sub); }};
                
                
        tilde_list.push_back(matched_tilde);
        
        auto finish2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed2 = finish2 - start2;
        std::cout<<"Done in "<< elapsed2.count()  << " seconds!" <<'\n';
        
        //-- Basis Amplitude -------------------------------------------
        auto start3 = std::chrono::high_resolution_clock::now();
        std::cout<<"Initializing BasisAmplitude ...\n";
        
        // Initialize the correct basisfunction
        basis::BasisAmplitudeEta3Pi amp(mass_1,decay_mass,
                                          omn,
                                          tilde_list[j],
                                          phases, cutoff, 
                                          epsilon,validity,max_subs);
                
        const auto basis_amp{[amp](double s, int isospin, SubtractionConstant sub, Setting evaluation){
                    return amp.amplitude(s, isospin, sub, evaluation); }};
                
        basis_list.push_back(basis_amp);
        
        std::cout<<"Finsihed initializing BasisAmplitude ...\n";
        
        auto finish3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed3 = finish3 - start3;
        std::cout<<"Done in "<< elapsed3.count()  << " seconds!" <<'\n';

        std::chrono::duration<double> elapsedz = finish3 - start2;
        std::cout<<"Finished iteration "<< j+1<< " in a total of " << elapsedz.count()  << " seconds!" <<'\n';
    }

    std::cout<<
    "\n----------------------------\n"
    "Finished all iterations!"
    "\n----------------------------\n";
    
    std::cout<<"Writing lists for final solution ...\n";
    write_output(SubtractionConstant::a0, iterations, output_a0);
    write_output(SubtractionConstant::b0, iterations, output_b0);
    write_output(SubtractionConstant::a1, iterations, output_a1);
    write_output(SubtractionConstant::g1, iterations, output_g1);
    write_output(SubtractionConstant::h1, iterations, output_h1);

    auto finish4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed4 = finish4 - start1;
    std::cout<<
    "Solved Khuri-Treiman equations in a total of "<< elapsed4.count()  << " seconds!" <<
    "\n------------------------------------------------------------------------------------------\n"<<
    "\n------------------------------------------------------------------------------------------\n\n";
    return;
}

void IterationEtap3Pi::solve_khuri_treiman(){
    if ((max_subs != 0) ){
        throw std::domain_error{"Etap -> 3Pi is only implemented for a constant subtraction scheme. Set max_subs to zero."};
    }
    auto start1 = std::chrono::high_resolution_clock::now();

    const auto omn
    {[this](double s, int isospin, Setting evaluation){
        return evaluate_omnes(s,isospin,evaluation);
    }};    

    auto starthom = std::chrono::high_resolution_clock::now();
    std::cout<<"Building Homogeneous solution ...\n";

    basis::HomogeneousSolutionEtap3Pi homsol(mass_1,decay_mass, omn);
    
    const auto hom{[homsol](double s, int isospin, SubtractionConstant sub, Setting evaluation){
        return homsol.initial_guess(s, isospin, sub, evaluation);}};
    
    basis_list.push_back(hom);
    
    auto finishhom = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedhom = finishhom - starthom;
    std::cout<<"Finished in "<< elapsedhom.count()  << " seconds!" <<'\n';
    
    for (int j=0; j<iterations; j++) {
        std::cout<<"----------- Iteration "<< j+1 << " ------------"<<"\n";
        //-- Inhomogeneity ---------------------------------------------
        auto start2 = std::chrono::high_resolution_clock::now();
        std::cout<<"Initializing BasisTilde ...\n";
       
        // Initialize the correct basisfunction
        basis::BasisTildeEtap3Pi tild(mass_1,decay_mass,
                                        basis_list[j], cutoff, j,
                                        epsilon_tilde, validity_tilde);
                
        const auto matched_tilde{[tild](double s, int isospin, SubtractionConstant sub){
            return tild.matched_tilde(s, isospin, sub); }};
                
                
        tilde_list.push_back(matched_tilde);

        
        auto finish2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed2 = finish2 - start2;
        std::cout<<"Done in "<< elapsed2.count()  << " seconds!" <<'\n';
        
        //-- Basis Amplitude -------------------------------------------
        auto start3 = std::chrono::high_resolution_clock::now();
        std::cout<<"Initializing BasisAmplitude ...\n";
        
        // Initialize the correct basisfunction
        basis::BasisAmplitudeEtap3Pi amp(mass_1,decay_mass,
                                          omn,
                                          tilde_list[j],
                                          phases, cutoff, 
                                          epsilon,validity,max_subs);
                
        const auto basis_amp{[amp](double s, int isospin, SubtractionConstant sub, Setting evaluation){
                    return amp.amplitude(s, isospin, sub, evaluation); }};
                
        basis_list.push_back(basis_amp);
        
        std::cout<<"Finsihed initializing BasisAmplitude ...\n";
        
        auto finish3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed3 = finish3 - start3;
        std::cout<<"Done in "<< elapsed3.count()  << " seconds!" <<'\n';

        std::chrono::duration<double> elapsedz = finish3 - start2;
        std::cout<<"Finished iteration "<< j+1<< " in a total of " << elapsedz.count()  << " seconds!" <<'\n';
    }

    std::cout<<
    "\n----------------------------\n"
    "Finished all iterations!"
    "\n----------------------------\n";
    
    std::cout<<"Writing lists for final solution ...\n";
    write_output(SubtractionConstant::a0, iterations, output_a0);
    write_output(SubtractionConstant::b0, iterations, output_b0);
    write_output(SubtractionConstant::a1, iterations, output_a1);
    write_output(SubtractionConstant::g1, iterations, output_g1);
    write_output(SubtractionConstant::h1, iterations, output_h1);

    auto finish4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed4 = finish4 - start1;
    std::cout<<
    "Solved Khuri-Treiman equations in a total of "<< elapsed4.count()  << " seconds!" <<
    "\n------------------------------------------------------------------------------------------\n"<<
    "\n------------------------------------------------------------------------------------------\n\n";
    return;
}

void IterationEtapEtaPiPi::solve_khuri_treiman(){
    if ((max_subs != 0) ){
        throw std::domain_error{"Etap -> Eta Pi Pi is only implemented for a constant subtraction scheme. Set max_subs to zero."};
    }
    auto start1 = std::chrono::high_resolution_clock::now();

    const auto omn
    {[this](double s, int isospin, Setting evaluation){
        return evaluate_omnes(s,isospin,evaluation);
    }};    

    auto starthom = std::chrono::high_resolution_clock::now();
    std::cout<<"Building Homogeneous solution ...\n";

    basis::HomogeneousSolutionEtapEtaPiPi homsol(mass_1,mass_3, decay_mass, omn);
    
    const auto hom{[homsol](double s, int isospin, SubtractionConstant sub, Setting evaluation){
        return homsol.initial_guess(s, isospin, sub, evaluation);}};
    
    basis_list.push_back(hom);
    
    auto finishhom = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedhom = finishhom - starthom;
    std::cout<<"Finished in "<< elapsedhom.count()  << " seconds!" <<'\n';
    
    for (int j=0; j<iterations; j++) {
        std::cout<<"----------- Iteration "<< j+1 << " ------------"<<"\n";
        //-- Inhomogeneity ---------------------------------------------
        auto start2 = std::chrono::high_resolution_clock::now();
        std::cout<<"Initializing BasisTilde ...\n";
       
        basis::BasisTildeEtapEtaPiPi tild(mass_1,mass_3,decay_mass,
                                basis_list[j], cutoff, j, 
                                epsilon_tilde, validity_tilde);
        
        const auto matched_tilde{[tild](double s, int isospin, SubtractionConstant sub){
            return tild.matched_tilde(s, isospin, sub); }};
        
        tilde_list.push_back(matched_tilde);
        
        auto finish2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed2 = finish2 - start2;
        std::cout<<"Done in "<< elapsed2.count()  << " seconds!" <<'\n';
        
        //-- Basis Amplitude -------------------------------------------
        auto start3 = std::chrono::high_resolution_clock::now();
        std::cout<<"Initializing BasisAmplitude ...\n";
        
        // Initialize the correct basisfunction
        basis::BasisAmplitudeEtapEtaPiPi amp(mass_1,mass_3,decay_mass,
                                  omn,
                                  tilde_list[j], phases, cutoff,
                                  epsilon,validity,max_subs);
        
        const auto basis_amp{[amp](double s, int isospin, SubtractionConstant sub, Setting evaluation){
            return amp.amplitude(s, isospin, sub, evaluation); }};
        
        basis_list.push_back(basis_amp);        
        
        std::cout<<"Finsihed initializing BasisAmplitude ...\n";
        
        auto finish3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed3 = finish3 - start3;
        std::cout<<"Done in "<< elapsed3.count()  << " seconds!" <<'\n';

        std::chrono::duration<double> elapsedz = finish3 - start2;
        std::cout<<"Finished iteration "<< j+1<< " in a total of " << elapsedz.count()  << " seconds!" <<'\n';
    }

    std::cout<<
    "\n----------------------------\n"
    "Finished all iterations!"
    "\n----------------------------\n";
    
    std::cout<<"Writing lists for final solution ...\n";
    write_output(SubtractionConstant::sub_a0, iterations, output_sub_a0);
    write_output(SubtractionConstant::sub_b0, iterations, output_sub_b0);
    write_output(SubtractionConstant::sub_c0, iterations, output_sub_c0);
    write_output(SubtractionConstant::sub_d0, iterations, output_sub_d0);
    write_output(SubtractionConstant::sub_a1, iterations, output_sub_a1);
    write_output(SubtractionConstant::sub_b1, iterations, output_sub_b1);
    
    auto finish4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed4 = finish4 - start1;
    std::cout<<
    "Solved Khuri-Treiman equations in a total of "<< elapsed4.count()  << " seconds!" <<
    "\n------------------------------------------------------------------------------------------\n"<<
    "\n------------------------------------------------------------------------------------------\n\n";
    return;
}

void IterationV3Pi::solve_khuri_treiman(){
    if ((max_subs != 1) && (max_subs != 2)){
        throw std::domain_error{"V -> 3Pi is only implemented for one or two subtractions."};
    }

    auto start1 = std::chrono::high_resolution_clock::now();

    const auto omn
    {[this](double s, int isospin, Setting evaluation){
        return evaluate_omnes(s,isospin,evaluation);
    }};    

    // std::cout << "Omnes at s=0" << omn(0,1,Setting::above) << std::endl;

    auto starthom = std::chrono::high_resolution_clock::now();
    std::cout<<"Building Homogeneous solution ...\n";

    basis::HomogeneousSolutionV3Pi homsol(mass_1,decay_mass, omn);
    
    const auto hom{[homsol](double s, int isospin, SubtractionConstant sub, Setting evaluation){
        return homsol.initial_guess(s, isospin, sub, evaluation);}};
    
    basis_list.push_back(hom);
    
    auto finishhom = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedhom = finishhom - starthom;
    std::cout<<"Finished in "<< elapsedhom.count()  << " seconds!" <<'\n';
    
    for (int j=0; j<iterations; j++) {
        std::cout<<"----------- Iteration "<< j+1 << " ------------"<<"\n";
        //-- Inhomogeneity ---------------------------------------------
        auto start2 = std::chrono::high_resolution_clock::now();
        std::cout<<"Initializing BasisTilde ...\n";
       
        // Initialize the correct basisfunction
        basis::BasisTildeV3Pi tild(mass_1,decay_mass,
                                        basis_list[j], cutoff, j,
                                        epsilon_tilde, validity_tilde);
                
        const auto matched_tilde{[tild](double s, int isospin, SubtractionConstant sub){
            return tild.matched_tilde(s, isospin, sub); }};
                
                
        tilde_list.push_back(matched_tilde);

        
        auto finish2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed2 = finish2 - start2;
        std::cout<<"Done in "<< elapsed2.count()  << " seconds!" <<'\n';
        
        //-- Basis Amplitude -------------------------------------------
        auto start3 = std::chrono::high_resolution_clock::now();
        std::cout<<"Initializing BasisAmplitude ...\n";
        
        // Initialize the correct basisfunction
        basis::BasisAmplitudeV3Pi amp(mass_1,decay_mass,
                                          omn,
                                          tilde_list[j],
                                          phases, cutoff, 
                                          epsilon,validity,max_subs);
                
        const auto basis_amp{[amp](double s, int isospin, SubtractionConstant sub, Setting evaluation){
                    return amp.amplitude(s, isospin, sub, evaluation); }};
                
        basis_list.push_back(basis_amp);
        
        std::cout<<"Finsihed initializing BasisAmplitude ...\n";
        
        auto finish3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed3 = finish3 - start3;
        std::cout<<"Done in "<< elapsed3.count()  << " seconds!" <<'\n';

        std::chrono::duration<double> elapsedz = finish3 - start2;
        std::cout<<"Finished iteration "<< j+1<< " in a total of " << elapsedz.count()  << " seconds!" <<'\n';
    }

    std::cout<<
    "\n----------------------------\n"
    "Finished all iterations!"
    "\n----------------------------\n";
    
    std::cout<<"Writing lists for final solution ...\n";
    if (max_subs == 1) {
        write_output(SubtractionConstant::a1, iterations, output_a);
    }
    else{
        write_output(SubtractionConstant::a1, iterations, output_a);
        write_output(SubtractionConstant::b1, iterations, output_b);
    }

    auto finish4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed4 = finish4 - start1;
    std::cout<<
    "Solved Khuri-Treiman equations in a total of "<< elapsed4.count()  << " seconds!" <<
    "\n------------------------------------------------------------------------------------------\n"<<
    "\n------------------------------------------------------------------------------------------\n\n";
    return;
}

void IterationX3Pi::solve_khuri_treiman(){
    if ((max_subs != 1) && (max_subs != 2)){
        throw std::domain_error{"V -> 3Pi is only implemented for one or two subtractions."};
    }

    auto start1 = std::chrono::high_resolution_clock::now();

    const auto omn
    {[this](double s, int isospin, Setting evaluation){
        return evaluate_omnes(s,isospin,evaluation);
    }};    

    auto starthom = std::chrono::high_resolution_clock::now();
    std::cout<<"Building Homogeneous solution ...\n";

    basis::HomogeneousSolutionX3Pi homsol(mass_1,decay_mass, omn);
    
    const auto hom{[homsol](double s, int isospin, SubtractionConstant sub, Setting evaluation){
        return homsol.initial_guess(s, isospin, sub, evaluation);}};
    
    basis_list.push_back(hom);
    
    auto finishhom = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedhom = finishhom - starthom;
    std::cout<<"Finished in "<< elapsedhom.count()  << " seconds!" <<'\n';
    
    for (int j=0; j<iterations; j++) {
        std::cout<<"----------- Iteration "<< j+1 << " ------------"<<"\n";
        //-- Inhomogeneity ---------------------------------------------
        auto start2 = std::chrono::high_resolution_clock::now();
        std::cout<<"Initializing BasisTilde ...\n";
       
        // Initialize the correct basisfunction
        basis::BasisTildeX3Pi tild(mass_1,decay_mass,
                                        basis_list[j], cutoff, j,
                                        epsilon_tilde, validity_tilde);
                
        const auto matched_tilde{[tild](double s, int isospin, SubtractionConstant sub){
            return tild.matched_tilde(s, isospin, sub); }};
                
                
        tilde_list.push_back(matched_tilde);

        
        auto finish2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed2 = finish2 - start2;
        std::cout<<"Done in "<< elapsed2.count()  << " seconds!" <<'\n';
        
        //-- Basis Amplitude -------------------------------------------
        auto start3 = std::chrono::high_resolution_clock::now();
        std::cout<<"Initializing BasisAmplitude ...\n";
        
        // Initialize the correct basisfunction
        basis::BasisAmplitudeX3Pi amp(mass_1,decay_mass,
                                          omn,
                                          tilde_list[j],
                                          phases, cutoff, 
                                          epsilon,validity,max_subs);
                
        const auto basis_amp{[amp](double s, int isospin, SubtractionConstant sub, Setting evaluation){
                    return amp.amplitude(s, isospin, sub, evaluation); }};
                
        basis_list.push_back(basis_amp);
        
        std::cout<<"Finsihed initializing BasisAmplitude ...\n";
        
        auto finish3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed3 = finish3 - start3;
        std::cout<<"Done in "<< elapsed3.count()  << " seconds!" <<'\n';

        std::chrono::duration<double> elapsedz = finish3 - start2;
        std::cout<<"Finished iteration "<< j+1<< " in a total of " << elapsedz.count()  << " seconds!" <<'\n';
    }

    std::cout<<
    "\n----------------------------\n"
    "Finished all iterations!"
    "\n----------------------------\n";
    
    std::cout<<"Writing lists for final solution ...\n";
    if (max_subs == 1) {
        write_output(SubtractionConstant::a1, iterations, output_a);
    }
    else{
        write_output(SubtractionConstant::a1, iterations, output_a);
        write_output(SubtractionConstant::b1, iterations, output_b);
    }

    auto finish4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed4 = finish4 - start1;
    std::cout<<
    "Solved Khuri-Treiman equations in a total of "<< elapsed4.count()  << " seconds!" <<
    "\n------------------------------------------------------------------------------------------\n"<<
    "\n------------------------------------------------------------------------------------------\n\n";
    return;
}

void IterationT3Pi::solve_khuri_treiman(){
    if ((max_subs != 1) && (max_subs != 2)){
        throw std::domain_error{"T -> 3Pi is only implemented for one or two subtractions."};
    }

    auto start1 = std::chrono::high_resolution_clock::now();

    const auto omn
    {[this](double s, int isospin, Setting evaluation){
        return evaluate_omnes(s,isospin,evaluation);
    }};    

    auto starthom = std::chrono::high_resolution_clock::now();
    std::cout<<"Building Homogeneous solution ...\n";

    basis::HomogeneousSolutionT3Pi homsol(mass_1,decay_mass, omn);
    
    const auto hom{[homsol](double s, int isospin, SubtractionConstant sub, Setting evaluation){
        return homsol.initial_guess(s, isospin, sub, evaluation);}};
    
    basis_list.push_back(hom);
    
    auto finishhom = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsedhom = finishhom - starthom;
    std::cout<<"Finished in "<< elapsedhom.count()  << " seconds!" <<'\n';
    
    for (int j=0; j<iterations; j++) {
        std::cout<<"----------- Iteration "<< j+1 << " ------------"<<"\n";
        //-- Inhomogeneity ---------------------------------------------
        auto start2 = std::chrono::high_resolution_clock::now();
        std::cout<<"Initializing BasisTilde ...\n";
       
        // Initialize the correct basisfunction
        basis::BasisTildeT3Pi tild(mass_1,decay_mass,
                                        basis_list[j], cutoff, j,
                                        epsilon_tilde, validity_tilde);
                
        const auto matched_tilde{[tild](double s, int isospin, SubtractionConstant sub){
            return tild.matched_tilde(s, isospin, sub); }};
                
                
        tilde_list.push_back(matched_tilde);

        
        auto finish2 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed2 = finish2 - start2;
        std::cout<<"Done in "<< elapsed2.count()  << " seconds!" <<'\n';
        
        //-- Basis Amplitude -------------------------------------------
        auto start3 = std::chrono::high_resolution_clock::now();
        std::cout<<"Initializing BasisAmplitude ...\n";
        
        // Initialize the correct basisfunction
        basis::BasisAmplitudeT3Pi amp(mass_1,decay_mass,
                                          omn,
                                          tilde_list[j],
                                          phases, cutoff, 
                                          epsilon,validity,max_subs);
                
        const auto basis_amp{[amp](double s, int isospin, SubtractionConstant sub, Setting evaluation){
                    return amp.amplitude(s, isospin, sub, evaluation); }};
                
        basis_list.push_back(basis_amp);
        
        std::cout<<"Finsihed initializing BasisAmplitude ...\n";
        
        auto finish3 = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> elapsed3 = finish3 - start3;
        std::cout<<"Done in "<< elapsed3.count()  << " seconds!" <<'\n';

        std::chrono::duration<double> elapsedz = finish3 - start2;
        std::cout<<"Finished iteration "<< j+1<< " in a total of " << elapsedz.count()  << " seconds!" <<'\n';
    }

    std::cout<<
    "\n----------------------------\n"
    "Finished all iterations!"
    "\n----------------------------\n";
    
    std::cout<<"Writing lists for final solution ...\n";
    if (max_subs == 1) {
        write_output(SubtractionConstant::a1, iterations, output_a);
    }
    else{
        write_output(SubtractionConstant::a1, iterations, output_a);
        write_output(SubtractionConstant::b1, iterations, output_b);
    }

    auto finish4 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed4 = finish4 - start1;
    std::cout<<
    "Solved Khuri-Treiman equations in a total of "<< elapsed4.count()  << " seconds!" <<
    "\n------------------------------------------------------------------------------------------\n"<<
    "\n------------------------------------------------------------------------------------------\n\n";
    return;
}


Complex Iteration::operator()(double s, int isospin, SubtractionConstant sub, Setting set, int iteration){
    return basis_list[iteration](s,isospin,sub,set);
}

}// namespace
