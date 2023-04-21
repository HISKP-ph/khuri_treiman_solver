#ifndef ITERATIVE_SOLUTION_H
#define ITERATIVE_SOLUTION_H

#include "phase.h"
#include "omnes.h"
#include "path.h"
#include "cauchy.h"
#include "gsl_interface.h"
#include "angular_average.h"
#include "enums.h"
#include "basis_function.h"
#include "type_aliases.h"
#include "constants.h"
#include "splined_omnes.h"
#include "splined_path.h"


#include <vector>
#include <string>
#include <fstream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

#include<complex>
#include<cmath>
#include<iostream>
#include<iomanip>
#include<chrono>

using type_aliases::Complex;
using type_aliases::CFunction;
using namespace enums;
using namespace std::complex_literals;

///Solve the Khuri-Treiman equations iteratively. 
namespace iteration{

///Parent class for the iterative solution of the Khuri-Treiman equations.
class Iteration{

public:
    Iteration(double mass_1, double mass_2,
              double mass_3, double mass_4,
              double epsilon, double validity,
              double epsilon_tilde, double validity_tilde,
              double cutoff, int iterations, int max_subs);
    ///<@param mass_1, mass_2, mass_3 masses of the decay products
    ///<@param mass_4 decay mass
    ///<@param epsilon, validity matching point and range of the expansion (around the pseudo-threshold) within the dispersion integral is used.
    ///<@param epsilon_tilde, validity_tilde matching point and range of the expansion (around the remaining singularities) within the dispersion integral is used.
    ///<@param cutoff numeric cutoff of the dispersion integral
    ///<@param iterations number of iterations
    ///<@param max_subs maximum number of subtraction constants (not valid for all of the processes)
    
    double mass_1, mass_2, mass_3, decay_mass;
    double epsilon, validity, epsilon_tilde, validity_tilde, cutoff;
    int iterations;
    int max_subs;
    /// This list contains the result of the angular average for each iterative step.
    std::vector<std::function<Complex (double, int, SubtractionConstant)>> tilde_list;
    /// This list contains the result for the basis amplitudes for each iterative step.
    std::vector<std::function<Complex (double, int, SubtractionConstant, Setting)>> basis_list;

    /// Function that solves the KT equations. Includes timings for the individual steps.
    virtual void solve_khuri_treiman()=0;
    /// Function that returns the Omnès function. (To easily use this in the python module and compare to the KT.)
    virtual Complex evaluate_omnes(double s, int isospin, Setting evaluation)=0;
    /// Returns the basis function of at the desired iteration.
    Complex operator()(double s, int isospin, SubtractionConstant sub, Setting set, int iteration);    
    /// Writes solution after this iteration into a file with name output_file. 
    virtual void write_output(SubtractionConstant sub, int iteration, const std::string output_file)=0;
    
};

class IterationEta3Pi: public Iteration{        
public:
    IterationEta3Pi(double mass_1, double mass_4,
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
              const std::string& output_h1);
    

    void solve_khuri_treiman() override;
    Complex evaluate_omnes(double s, int isospin, Setting evaluation) override;

private:   
    const std::string output_a0;
    const std::string output_b0;
    const std::string output_a1;
    const std::string output_g1;
    const std::string output_h1;
public:
    array::ArrayEta3Pi mandelstam_array;
private:
    phase::PhaseEta3Pi phases;
    splined_omnes::OmnesSplineEta3Pi omnes;
    
    void write_output(SubtractionConstant, int iteration, const std::string output_file) override;
};

class IterationEtap3Pi: public Iteration{        
public:
    IterationEtap3Pi(double mass_1, double mass_4,
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
              const std::string& output_h1);
    

    void solve_khuri_treiman() override;
    Complex evaluate_omnes(double s, int isospin, Setting evaluation) override;
private:    
    const std::string output_a0;
    const std::string output_b0;
    const std::string output_a1;
    const std::string output_g1;
    const std::string output_h1;
    
    array::ArrayEtap3Pi mandelstam_array;
    phase::PhaseEtap3Pi phases;
    splined_omnes::OmnesSplineEtap3Pi omnes;
    
    void write_output(SubtractionConstant sub, int iteration, const std::string output_file) override;
};

class IterationEtapEtaPiPi:public Iteration{
public:
    IterationEtapEtaPiPi(double mass_1,
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
              const std::string& output_sub_b1);
    

    void solve_khuri_treiman() override;
    Complex evaluate_omnes(double s, int isospin, Setting evaluation) override;
private:   
    const std::string output_sub_a0;
    const std::string output_sub_b0;
    const std::string output_sub_c0;
    const std::string output_sub_d0;
    const std::string output_sub_a1;
    const std::string output_sub_b1;
    
    array::ArrayEtapEtaPiPi mandelstam_array;
    phase::PhaseEtapEtaPiPi phases;
    splined_omnes::OmnesSplineEtapEtaPiPi omnes;
    
    void write_output(SubtractionConstant sub, int iteration, const std::string output_file) override;
};

class IterationV3Pi: public Iteration{        
public:
    IterationV3Pi(double mass_1, double mass_4,
              const std::string& phase,
              double epsilon, double validity,
              double epsilon_tilde, double validity_tilde,
              double cutoff, int iterations, int max_subs,
              const std::string& output_a,
              const std::string& output_b);


    void solve_khuri_treiman() override;
    Complex evaluate_omnes(double s, int isospin, Setting evaluation) override;
private:    
    const std::string output_a;
    const std::string output_b;

public:
    array::ArrayV3Pi mandelstam_array;  
private:
    phase::PhaseV3Pi phases;
    splined_omnes::OmnesSplineV3Pi omnes;
    
    void write_output(SubtractionConstant sub, int iteration, const std::string output_file) override;
};

class IterationX3Pi: public Iteration{        
public:
    IterationX3Pi(double mass_1, double mass_4,
              const std::string& phase,
              double epsilon, double validity,
              double epsilon_tilde, double validity_tilde,
              double cutoff, int iterations, int max_subs,
              const std::string& output_a,
              const std::string& output_b);


    void solve_khuri_treiman() override;
    Complex evaluate_omnes(double s, int isospin, Setting evaluation) override;
private:    
    const std::string output_a;
    const std::string output_b;

public:
    array::ArrayX3Pi mandelstam_array;  
private:
    phase::PhaseX3Pi phases;
    splined_omnes::OmnesSplineX3Pi omnes;
    
    void write_output(SubtractionConstant sub, int iteration, const std::string output_file) override;
};

class IterationT3Pi: public Iteration{        
public:
    IterationT3Pi(double mass_1, double mass_4,
              const std::string& phase,
              double epsilon, double validity,
              double epsilon_tilde, double validity_tilde,
              double cutoff, int iterations, int max_subs,
              const std::string& output_a,
              const std::string& output_b);


    void solve_khuri_treiman() override;
    Complex evaluate_omnes(double s, int isospin, Setting evaluation) override;
private:    
    const std::string output_a;
    const std::string output_b;

public:
    array::ArrayT3Pi mandelstam_array;  
private:
    phase::PhaseT3Pi phases;
    splined_omnes::OmnesSplineT3Pi omnes;
    
    void write_output(SubtractionConstant sub, int iteration, const std::string output_file) override;
};


} // iteration

#endif // ITERATIVE_SOLUTION_H
