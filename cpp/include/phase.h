// Read in the phase shifts from an arbitrary file and interpolate the respective lists.
// Attention: the lists have to use Mandelstam s in units of the squared pion mass

#ifndef PHASE_H
#define PHASE_H

#include <vector>
#include <cmath>
#include <fstream>
#include <string>
#include "gsl_interface.h"
#include "asymptotic.h"


using gsl::Interpolate;

namespace phase{

/// Parent class for reading, splining and continuing the phase shift.
class Phase{
public:
    Phase() = default;
    /// Returns the phase at mandelstam s for a defined isospin.
    virtual double phase(double s, int isospin)const=0;
    /// Reads in a file with name phase containing a list of mandelstam s and the phase shift.  
    void readin_phase(const std::string& phase,
        std::vector<double>& s_list, std::vector<double>& delta_list);
private:
    /// Splines the read in phase.
    virtual void spline_phases()=0;
    /// Continues the phase to higher energies.
    virtual void make_cont()=0;
};

class PhaseEta3Pi: public Phase{
public:
    PhaseEta3Pi(std::string phase_0, std::string phase_1, std::string phase_2);
    ///<@param phase_I the complete path name of the phase shift with isospin I that shall be imported
    
    double phase(double s, int isospin)const override;
    
private:
    void spline_phases() override;
    void make_cont() override;

    /// list of Mandelstam s for I=0
    std::vector<double> s_list_0;
    /// list of pi-pi  phase shift for I=0
    std::vector<double> delta_list_0; 
    
    /// list of Mandelstam s for I=1
    std::vector<double> s_list_1;
    /// list of pi-pi phase shift for I=1
    std::vector<double> delta_list_1; 
    
    /// list of Mandelstam s for I=2
    std::vector<double> s_list_2;
    /// list of pi-pi phase shift for I=2
    std::vector<double> delta_list_2; 
    
    Interpolate spline_phase_0;
    Interpolate spline_phase_1;
    Interpolate spline_phase_2;

    asymptotic::Asymptotic1s cont0;
    asymptotic::Asymptotic1s cont1;   
};

class PhaseEtap3Pi: public Phase{
public:
    PhaseEtap3Pi(std::string phase_0, std::string phase_1, std::string phase_2);
    ///<@param phase_I the complete path name of the phase shift with isospin I that shall be imported
    
    double phase(double s, int isospin)const override;
    
private:
    void spline_phases() override;
    void make_cont() override;

    /// list of Mandelstam s for I=0
    std::vector<double> s_list_0;
    /// list of pi-pi  phase shift for I=0
    std::vector<double> delta_list_0; 
    
    /// list of Mandelstam s for I=1
    std::vector<double> s_list_1;
    /// list of pi-pi phase shift for I=1
    std::vector<double> delta_list_1; 
    
    /// list of Mandelstam s for I=2
    std::vector<double> s_list_2;
    /// list of pi-pi phase shift for I=2
    std::vector<double> delta_list_2; 
    
    Interpolate spline_phase_0;
    Interpolate spline_phase_1;
    Interpolate spline_phase_2;

    asymptotic::Asymptotic1s cont0;
    asymptotic::Asymptotic1s cont1;  
};

class PhaseEtapEtaPiPi: public Phase{
public:
    PhaseEtapEtaPiPi(std::string phase_0, std::string phase_1, std::string phase_2, std::string phase_eta_pi);
    ///<@param phase_I the complete path name of the phase shift with isospin I that shall be imported
    
    double phase(double s, int isospin)const override;
    
    /// S-wave eta-pi phase shift
    double phase_eta_pi(double s)const;
    
private:
    void spline_phases() override;
    void make_cont() override;

    /// list of Mandelstam s for I=0
    std::vector<double> s_list_0;
    /// list of pi-pi  phase shift for I=0
    std::vector<double> delta_list_0; 
    
    /// list of Mandelstam s for I=1
    std::vector<double> s_list_1;
    /// list of pi-pi phase shift for I=1
    std::vector<double> delta_list_1; 
    
    /// list of Mandelstam s for I=2
    std::vector<double> s_list_2;
    /// list of pi-pi phase shift for I=2
    std::vector<double> delta_list_2; 
    
    /// list of Mandelstam s for eta-pi phase shift
    std::vector<double> s_list_eta_pi;
    /// list of for eta-pi  phase shift   
    std::vector<double> delta_list_eta_pi; 
    
    Interpolate spline_phase_0;
    Interpolate spline_phase_1;
    Interpolate spline_phase_2;
    Interpolate spline_phase_eta_pi;

    asymptotic::Asymptotic1s cont0;
    asymptotic::Asymptotic1s cont1;
    asymptotic::Asymptotic1s cont_eta_pi;  
};

class PhaseV3Pi: public Phase{
public:
    PhaseV3Pi(std::string phase);
    ///<@param phase_I the complete path name of the phase shift with isospin I that shall be imported
    
    double phase(double s, int isospin)const override;
private:
    void spline_phases() override;
    void make_cont() override;

    /// list of Mandelstam s
    std::vector<double> s_list;
    /// list of pi-pi phase shift
    std::vector<double> delta_list; 
    
    Interpolate spline_phase;

    asymptotic::Asymptotic1s cont;
};

class PhaseX3Pi: public Phase{
public:
    PhaseX3Pi(std::string phase);
    ///<@param phase_I the complete path name of the phase shift with isospin I that shall be imported
    
    double phase(double s, int isospin)const override;

private:
    void spline_phases() override;
    void make_cont() override;

    /// list of Mandelstam s
    std::vector<double> s_list;
    /// list of pi-pi phase shift
    std::vector<double> delta_list; 
    
    Interpolate spline_phase;

    asymptotic::Asymptotic1s cont;
};

class PhaseT3Pi: public Phase{
public:
    PhaseT3Pi(std::string phase);
    ///<@param phase_I the complete path name of the phase shift with isospin I that shall be imported
    
    double phase(double s, int isospin)const override;
    
private:
    void spline_phases() override;
    void make_cont() override;

    /// list of Mandelstam s
    std::vector<double> s_list;
    /// list of pi-pi phase shift
    std::vector<double> delta_list;  
    
    Interpolate spline_phase;

    asymptotic::Asymptotic1s cont;
};

} // phase

#endif // PHASE_H
