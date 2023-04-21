#include "phase.h"
#include "constants.h"
using constants::pi;

using namespace phase;



//-- Constructors -----------------------------------------------
PhaseEta3Pi::PhaseEta3Pi(std::string phase_0, std::string phase_1, std::string phase_2)
{readin_phase(phase_0,s_list_0,delta_list_0); readin_phase(phase_1,s_list_1,delta_list_1);
readin_phase(phase_2,s_list_2,delta_list_2); spline_phases(); make_cont();}

PhaseEtap3Pi::PhaseEtap3Pi(std::string phase_0, std::string phase_1, std::string phase_2)
{readin_phase(phase_0,s_list_0,delta_list_0); readin_phase(phase_1,s_list_1,delta_list_1);
readin_phase(phase_2,s_list_2,delta_list_2); spline_phases(); make_cont();}

PhaseEtapEtaPiPi::PhaseEtapEtaPiPi(std::string phase_0, std::string phase_1, std::string phase_2, std::string phase_eta_pi)
{readin_phase(phase_0,s_list_0,delta_list_0); readin_phase(phase_1,s_list_1,delta_list_1);
readin_phase(phase_2,s_list_2,delta_list_2); readin_phase(phase_eta_pi,s_list_eta_pi,delta_list_eta_pi); spline_phases(); make_cont();}

PhaseV3Pi::PhaseV3Pi(std::string phase)
{readin_phase(phase,s_list,delta_list); spline_phases(); make_cont();}

PhaseX3Pi::PhaseX3Pi(std::string phase)
{readin_phase(phase,s_list,delta_list); spline_phases(); make_cont();}

PhaseT3Pi::PhaseT3Pi(std::string phase)
{readin_phase(phase,s_list,delta_list); spline_phases(); make_cont();}
//--------------------------------------------------------------

void Phase::readin_phase(const std::string& phase,
                              std::vector<double>& s_list, std::vector<double>& delta_list){
    double s, delta;
    std::ifstream phase_shift{phase};
    while (phase_shift >> s >> delta) {
        s_list.push_back(s);
        delta_list.push_back(delta);
    }
    return;
}


void PhaseEta3Pi::spline_phases(){
    spline_phase_0 = Interpolate(s_list_0,delta_list_0,gsl::InterpolationMethod::cubic);
    spline_phase_1 = Interpolate(s_list_1,delta_list_1,gsl::InterpolationMethod::cubic);
    spline_phase_2 = Interpolate(s_list_2,delta_list_2,gsl::InterpolationMethod::cubic);
    return;
}

void PhaseEtap3Pi::spline_phases(){
    spline_phase_0 = Interpolate(s_list_0,delta_list_0,gsl::InterpolationMethod::cubic);
    spline_phase_1 = Interpolate(s_list_1,delta_list_1,gsl::InterpolationMethod::cubic);
    spline_phase_2 = Interpolate(s_list_2,delta_list_2,gsl::InterpolationMethod::cubic);
    return;
}

void PhaseEtapEtaPiPi::spline_phases(){
    spline_phase_0 = Interpolate(s_list_0,delta_list_0,gsl::InterpolationMethod::cubic);
    spline_phase_1 = Interpolate(s_list_1,delta_list_1,gsl::InterpolationMethod::cubic);
    spline_phase_2 = Interpolate(s_list_2,delta_list_2,gsl::InterpolationMethod::cubic);
    spline_phase_eta_pi = Interpolate(s_list_eta_pi,delta_list_eta_pi,gsl::InterpolationMethod::cubic);
    return;
}

void PhaseV3Pi::spline_phases(){
    spline_phase = Interpolate(s_list,delta_list,gsl::InterpolationMethod::cubic);
    return;
}

void PhaseX3Pi::spline_phases(){
    spline_phase = Interpolate(s_list,delta_list,gsl::InterpolationMethod::cubic);
    return;
}

void PhaseT3Pi::spline_phases(){
    spline_phase = Interpolate(s_list,delta_list,gsl::InterpolationMethod::cubic);
    return;
}

void PhaseEta3Pi::make_cont(){
    cont0 = asymptotic::Asymptotic1s(spline_phase_0, 115);
    cont1 = asymptotic::Asymptotic1s(spline_phase_1, 80);
}

void PhaseEtap3Pi::make_cont(){
    cont0 = asymptotic::Asymptotic1s(spline_phase_0, 115);
    cont1 = asymptotic::Asymptotic1s(spline_phase_1, 80);
}

void PhaseEtapEtaPiPi::make_cont(){
    cont0 = asymptotic::Asymptotic1s(spline_phase_0, 115);
    cont1 = asymptotic::Asymptotic1s(spline_phase_1, 80);
    cont_eta_pi = asymptotic::Asymptotic1s(spline_phase_eta_pi, 110);
}

void PhaseV3Pi::make_cont(){
    cont = asymptotic::Asymptotic1s(spline_phase, 80);
}

void PhaseX3Pi::make_cont(){
    cont = asymptotic::Asymptotic1s(spline_phase, 80);
}

void PhaseT3Pi::make_cont(){
    cont = asymptotic::Asymptotic1s(spline_phase, 80);
}

// Use an arbitrary continuation of the phase shifts for large energies.
// Here the continuation is only valid for the phase shifts of the
// 'Bern' analysis.
double PhaseEta3Pi::phase(double s, int isospin)const{
    if (isospin==0) {
        if (s<=115) {
            return spline_phase_0(s);
        }
        else if (s>115) {
            return cont0(s);
        }
        throw std::domain_error{"Value of mandelstam s not allowed."};
    }
    if (isospin==1) {
        if (s<=80) {
            return spline_phase_1(s);
        }
        else if (s>80) {
            return cont1(s);
        }
        throw std::domain_error{"Value of mandelstam s not allowed."};
    }
    if (isospin==2) {
        if (s<=800) {
            return spline_phase_2(s);
        }
        else if (s>800) {
            return 0;
        }
        throw std::domain_error{"Value of mandelstam s not allowed."};
    }
    else{throw std::domain_error{"The isospin entering 'phase' has to be either 0,1 or 2."};}
}


double PhaseEtap3Pi::phase(double s, int isospin)const{
    if (isospin==0) {
        if (s<=115) {
            return spline_phase_0(s);
        }
        else if (s>115) {
            return cont0(s);
        }
        throw std::domain_error{"Value of mandelstam s not allowed."};
    }
    if (isospin==1) {
        if (s<=80) {
            return spline_phase_1(s);
        }
        else if (s>80) {
            return cont1(s);
        }
        throw std::domain_error{"Value of mandelstam s not allowed."};
    }
    if (isospin==2) {
        if (s<=800) {
            return spline_phase_2(s);
        }
        else if (s>800) {
            return 0;
        }
        throw std::domain_error{"Value of mandelstam s not allowed."};
    }
    else{throw std::domain_error{"The isospin entering 'phase' has to be either 0,1 or 2."};}
}


double PhaseEtapEtaPiPi::phase(double s, int isospin)const{
    if (isospin==0) {
        if (s<=115) {
            return spline_phase_0(s);
        }
        else if (s>115) {
            return cont0(s);
        }
        throw std::domain_error{"Value of mandelstam s not allowed."};
    }
    if (isospin==1) {
        if (s<=80) {
            return spline_phase_1(s);
        }
        else if (s>80) {
            return cont1(s);
        }
        throw std::domain_error{"Value of mandelstam s not allowed."};
    }
    if (isospin==2) {
        if (s<=800) {
            return spline_phase_2(s);
        }
        else if (s>800) {
            return 0;
        }
        throw std::domain_error{"Value of mandelstam s not allowed."};
    }
    else{throw std::domain_error{"The isospin entering 'phase' has to be either 0,1 or 2."};}
}

double PhaseEtapEtaPiPi::phase_eta_pi(double s)const{
    if (s<=110) {
        return spline_phase_eta_pi(s);
    }
    else if (s>110) {
        return cont_eta_pi(s);
    }
    throw std::domain_error{"Value of mandelstam s not allowed."};
}

double PhaseV3Pi::phase(double s, int isospin)const{
    if (isospin==1) {
        if (s<=80) {
            return spline_phase(s);
        }
        else if (s>80) {
            return cont(s);
        }
        throw std::domain_error{"Value of mandelstam s not allowed."};
    }
    else{throw std::domain_error{"The isospin entering 'phase' has to be either 0,1 or 2."};}
}

double PhaseX3Pi::phase(double s, int isospin)const{
    if (isospin==1) {
        if (s<=80) {
            return spline_phase(s);
        }
        else if (s>80) {
            return cont(s);
        }
        throw std::domain_error{"Value of mandelstam s not allowed."};
    }
    else{throw std::domain_error{"The isospin entering 'phase' has to be either 0,1 or 2."};}
}

double PhaseT3Pi::phase(double s, int isospin)const{
    if (isospin==1) {
        if (s<=80) {
            return spline_phase(s);
        }
        else if (s>80) {
            return cont(s);
        }
        throw std::domain_error{"Value of mandelstam s not allowed."};
    }
    else{throw std::domain_error{"The isospin entering 'phase' has to be either 0,1 or 2."};}
}