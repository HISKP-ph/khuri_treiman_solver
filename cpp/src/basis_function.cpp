#include "basis_function.h"


namespace basis {


HomogeneousSolution::HomogeneousSolution(double mass_1, double mass_2, double mass_3, double mass_4,
                                         FunctionSetOmnes omnes)
:
omnes{omnes}
{}

HomogeneousSolutionEta3Pi::HomogeneousSolutionEta3Pi(double mass_1, double mass_4,
                                         FunctionSetOmnes omnes)
:
HomogeneousSolution(mass_1,mass_1,mass_1,mass_4,omnes),
polar_egg(mass_4, mass_1)
{}

HomogeneousSolutionEtap3Pi::HomogeneousSolutionEtap3Pi(double mass_1, double mass_4,
                                         FunctionSetOmnes omnes)
:
HomogeneousSolution(mass_1,mass_1,mass_1,mass_4,omnes),
polar_egg(mass_4, mass_1)
{}

HomogeneousSolutionEtapEtaPiPi::HomogeneousSolutionEtapEtaPiPi(double mass_1, double mass_3, double mass_4,
                                         FunctionSetOmnes omnes)
:
HomogeneousSolution(mass_1,mass_1,mass_3,mass_4,omnes),
gamma(mass_1,mass_1,mass_3,mass_4)
{}

HomogeneousSolutionV3Pi::HomogeneousSolutionV3Pi(double mass_1, double mass_4,
                                         FunctionSetOmnes omnes)
:
HomogeneousSolution(mass_1,mass_1,mass_1,mass_4,omnes),
polar_egg(mass_4, mass_1)
{}

HomogeneousSolutionX3Pi::HomogeneousSolutionX3Pi(double mass_1, double mass_4,
                                         FunctionSetOmnes omnes)
:
HomogeneousSolution(mass_1,mass_1,mass_1,mass_4,omnes),
polar_egg(mass_4, mass_1)
{}

HomogeneousSolutionT3Pi::HomogeneousSolutionT3Pi(double mass_1, double mass_4,
                                         FunctionSetOmnes omnes)
:
HomogeneousSolution(mass_1,mass_1,mass_1,mass_4,omnes),
polar_egg(mass_4, mass_1)
{}

BasisTildeEta3Pi::BasisTildeEta3Pi(double mass_1, double mass_4,
                       LongFunctionSet amplitude,
                       double cutoff, int iteration, double epsilon, double validity)
:
BasisTilde(),
// initialize angular anverages for each set of basis amplitudes and all possible isospins
ang_a0_0(mass_1,mass_4,amplitude, 0, cutoff, iteration, epsilon, validity, SubtractionConstant::a0),
ang_a0_1(mass_1,mass_4,amplitude, 1, cutoff, iteration, epsilon, validity, SubtractionConstant::a0),
ang_a0_2(mass_1,mass_4,amplitude, 2, cutoff, iteration, epsilon, validity, SubtractionConstant::a0),
ang_b0_0(mass_1,mass_4,amplitude, 0, cutoff, iteration, epsilon, validity, SubtractionConstant::b0),
ang_b0_1(mass_1,mass_4,amplitude, 1, cutoff, iteration, epsilon, validity, SubtractionConstant::b0),
ang_b0_2(mass_1,mass_4,amplitude, 2, cutoff, iteration, epsilon, validity, SubtractionConstant::b0),
ang_a1_0(mass_1,mass_4,amplitude, 0, cutoff, iteration, epsilon, validity, SubtractionConstant::a1),
ang_a1_1(mass_1,mass_4,amplitude, 1, cutoff, iteration, epsilon, validity, SubtractionConstant::a1),
ang_a1_2(mass_1,mass_4,amplitude, 2, cutoff, iteration, epsilon, validity, SubtractionConstant::a1),
//C-violating
ang_g1_1(mass_1,mass_4,amplitude, 1, cutoff, iteration, epsilon, validity, SubtractionConstant::g1),
ang_h1_1(mass_1,mass_4,amplitude, 1, cutoff, iteration, epsilon, validity, SubtractionConstant::h1),
ang_h1_2(mass_1,mass_4,amplitude, 2, cutoff, iteration, epsilon, validity, SubtractionConstant::h1)
{}

BasisTildeEtap3Pi::BasisTildeEtap3Pi(double mass_1, double mass_4,
                       LongFunctionSet amplitude,
                       double cutoff, int iteration, double epsilon, double validity)
:
BasisTilde(),
// initialize angular anverages for each set of basis amplitudes and all possible isospins
ang_a0_0(mass_1,mass_4,amplitude, 0, cutoff, iteration, epsilon, validity, SubtractionConstant::a0),
ang_a0_1(mass_1,mass_4,amplitude, 1, cutoff, iteration, epsilon, validity, SubtractionConstant::a0),
ang_a0_2(mass_1,mass_4,amplitude, 2, cutoff, iteration, epsilon, validity, SubtractionConstant::a0),
ang_b0_0(mass_1,mass_4,amplitude, 0, cutoff, iteration, epsilon, validity, SubtractionConstant::b0),
ang_b0_1(mass_1,mass_4,amplitude, 1, cutoff, iteration, epsilon, validity, SubtractionConstant::b0),
ang_b0_2(mass_1,mass_4,amplitude, 2, cutoff, iteration, epsilon, validity, SubtractionConstant::b0),
ang_a1_0(mass_1,mass_4,amplitude, 0, cutoff, iteration, epsilon, validity, SubtractionConstant::a1),
ang_a1_1(mass_1,mass_4,amplitude, 1, cutoff, iteration, epsilon, validity, SubtractionConstant::a1),
ang_a1_2(mass_1,mass_4,amplitude, 2, cutoff, iteration, epsilon, validity, SubtractionConstant::a1),
//C-violating
ang_g1_1(mass_1,mass_4,amplitude, 1, cutoff, iteration, epsilon, validity, SubtractionConstant::g1),
ang_h1_1(mass_1,mass_4,amplitude, 1, cutoff, iteration, epsilon, validity, SubtractionConstant::h1),
ang_h1_2(mass_1,mass_4,amplitude, 2, cutoff, iteration, epsilon, validity, SubtractionConstant::h1)
{}


BasisTildeEtapEtaPiPi::BasisTildeEtapEtaPiPi(double mass_1, double mass_3, double mass_4,
                       LongFunctionSet amplitude,
                       double cutoff, int iteration, double epsilon, double validity)
:
BasisTilde(),
// for eta'-> eta pi pi
/// Notation for the last three digits: e.g. ang_a0_XYZ has total isospin X, isospin of two-particle intermediate state Y and partial wave Z
ang_a0_000(mass_1,mass_3,mass_4,amplitude, 000, cutoff, iteration, epsilon, validity, SubtractionConstant::sub_a0),
ang_a0_010(mass_1,mass_3,mass_4,amplitude, 010, cutoff, iteration, epsilon, validity, SubtractionConstant::sub_a0),
ang_b0_000(mass_1,mass_3,mass_4,amplitude, 000, cutoff, iteration, epsilon, validity, SubtractionConstant::sub_b0),
ang_b0_010(mass_1,mass_3,mass_4,amplitude, 010, cutoff, iteration, epsilon, validity, SubtractionConstant::sub_b0),
ang_c0_000(mass_1,mass_3,mass_4,amplitude, 000, cutoff, iteration, epsilon, validity, SubtractionConstant::sub_c0),
ang_c0_010(mass_1,mass_3,mass_4,amplitude, 010, cutoff, iteration, epsilon, validity, SubtractionConstant::sub_c0),
ang_d0_000(mass_1,mass_3,mass_4,amplitude, 000, cutoff, iteration, epsilon, validity, SubtractionConstant::sub_d0),
ang_d0_010(mass_1,mass_3,mass_4,amplitude, 010, cutoff, iteration, epsilon, validity, SubtractionConstant::sub_d0),
//C-violating
ang_a1_111(mass_1,mass_3,mass_4,amplitude, 111, cutoff, iteration, epsilon, validity, SubtractionConstant::sub_a1),
ang_a1_110(mass_1,mass_3,mass_4,amplitude, 110, cutoff, iteration, epsilon, validity, SubtractionConstant::sub_a1),
ang_b1_111(mass_1,mass_3,mass_4,amplitude, 111, cutoff, iteration, epsilon, validity, SubtractionConstant::sub_b1),
ang_b1_110(mass_1,mass_3,mass_4,amplitude, 110, cutoff, iteration, epsilon, validity, SubtractionConstant::sub_b1)
{}

BasisTildeV3Pi::BasisTildeV3Pi(double mass_1, double mass_4,
                       LongFunctionSet amplitude,
                       double cutoff, int iteration, double epsilon, double validity)
:
BasisTilde(),
// initialize angular anverages for each set of basis amplitudes and all possible isospins
ang_a1(mass_1,mass_4,amplitude, 1, cutoff, iteration, epsilon, validity, SubtractionConstant::a1),
ang_b1(mass_1,mass_4,amplitude, 1, cutoff, iteration, epsilon, validity, SubtractionConstant::b1)
{}

BasisTildeX3Pi::BasisTildeX3Pi(double mass_1, double mass_4,
                       LongFunctionSet amplitude,
                       double cutoff, int iteration, double epsilon, double validity)
:
BasisTilde(),
// initialize angular anverages for each set of basis amplitudes and all possible isospins
ang_a1(mass_1,mass_4,amplitude, 1, cutoff, iteration, epsilon, validity, SubtractionConstant::a1),
ang_b1(mass_1,mass_4,amplitude, 1, cutoff, iteration, epsilon, validity, SubtractionConstant::b1)
{}

BasisTildeT3Pi::BasisTildeT3Pi(double mass_1, double mass_4,
                       LongFunctionSet amplitude,
                       double cutoff, int iteration, double epsilon, double validity)
:
BasisTilde(),
// initialize angular anverages for each set of basis amplitudes and all possible isospins
ang_a1(mass_1,mass_4,amplitude, 1, cutoff, iteration, epsilon, validity, SubtractionConstant::a1),
ang_b1(mass_1,mass_4,amplitude, 1, cutoff, iteration, epsilon, validity, SubtractionConstant::b1)
{}

BasisAmplitude::BasisAmplitude(double mass_1, double mass_2, double mass_3, double mass_4,
                               FunctionSetOmnes omnes,
                               FunctionSet matched_tilde,
                               double cutoff, double epsilon, double validity, int max_subs)
:
omnes{omnes},max_subs{max_subs}
{}

BasisAmplitudeEta3Pi::BasisAmplitudeEta3Pi(double mass_1, double mass_4,
                               FunctionSetOmnes omnes,
                               FunctionSet matched_tilde,
                               phase::PhaseEta3Pi phases,
                               double cutoff,double epsilon, double validity, int max_subs)
:
BasisAmplitude(mass_1,mass_1,mass_1,mass_4,omnes,matched_tilde,cutoff,epsilon,validity,max_subs),
polar_egg(mass_4, mass_1),
disp_a0_0(mass_1,mass_4,omnes,matched_tilde,phases, 0, epsilon, validity, cutoff, SubtractionConstant::a0, max_subs),
disp_a0_1(mass_1,mass_4,omnes,matched_tilde,phases, 1, epsilon, validity, cutoff, SubtractionConstant::a0, max_subs),
disp_a0_2(mass_1,mass_4,omnes,matched_tilde,phases, 2, epsilon, validity, cutoff, SubtractionConstant::a0, max_subs),
disp_b0_0(mass_1,mass_4,omnes,matched_tilde,phases, 0, epsilon, validity, cutoff, SubtractionConstant::b0, max_subs),
disp_b0_1(mass_1,mass_4,omnes,matched_tilde,phases, 1, epsilon, validity, cutoff, SubtractionConstant::b0, max_subs),
disp_b0_2(mass_1,mass_4,omnes,matched_tilde,phases, 2, epsilon, validity, cutoff, SubtractionConstant::b0, max_subs),
disp_a1_0(mass_1,mass_4,omnes,matched_tilde,phases, 0, epsilon, validity, cutoff, SubtractionConstant::a1, max_subs),
disp_a1_1(mass_1,mass_4,omnes,matched_tilde,phases, 1, epsilon, validity, cutoff, SubtractionConstant::a1, max_subs),
disp_a1_2(mass_1,mass_4,omnes,matched_tilde,phases, 2, epsilon, validity, cutoff, SubtractionConstant::a1, max_subs),
//C-violating
disp_g1_1(mass_1,mass_4,omnes,matched_tilde,phases, 1, epsilon, validity, cutoff, SubtractionConstant::g1, max_subs),
disp_h1_1(mass_1,mass_4,omnes,matched_tilde,phases, 1, epsilon, validity, cutoff, SubtractionConstant::h1, max_subs),
disp_h1_2(mass_1,mass_4,omnes,matched_tilde,phases, 2, epsilon, validity, cutoff, SubtractionConstant::h1, max_subs)
{}

BasisAmplitudeEtap3Pi::BasisAmplitudeEtap3Pi(double mass_1, double mass_4,
                               FunctionSetOmnes omnes,
                               FunctionSet matched_tilde,
                               phase::PhaseEtap3Pi phases,
                               double cutoff,double epsilon, double validity, int max_subs)
:
BasisAmplitude(mass_1,mass_1,mass_1,mass_4,omnes,matched_tilde,cutoff,epsilon,validity,max_subs),
polar_egg(mass_4, mass_1),
disp_a0_0(mass_1,mass_4,omnes,matched_tilde,phases, 0, epsilon, validity, cutoff, SubtractionConstant::a0, max_subs),
disp_a0_1(mass_1,mass_4,omnes,matched_tilde,phases, 1, epsilon, validity, cutoff, SubtractionConstant::a0, max_subs),
disp_a0_2(mass_1,mass_4,omnes,matched_tilde,phases, 2, epsilon, validity, cutoff, SubtractionConstant::a0, max_subs),
disp_b0_0(mass_1,mass_4,omnes,matched_tilde,phases, 0, epsilon, validity, cutoff, SubtractionConstant::b0, max_subs),
disp_b0_1(mass_1,mass_4,omnes,matched_tilde,phases, 1, epsilon, validity, cutoff, SubtractionConstant::b0, max_subs),
disp_b0_2(mass_1,mass_4,omnes,matched_tilde,phases, 2, epsilon, validity, cutoff, SubtractionConstant::b0, max_subs),
disp_a1_0(mass_1,mass_4,omnes,matched_tilde,phases, 0, epsilon, validity, cutoff, SubtractionConstant::a1, max_subs),
disp_a1_1(mass_1,mass_4,omnes,matched_tilde,phases, 1, epsilon, validity, cutoff, SubtractionConstant::a1, max_subs),
disp_a1_2(mass_1,mass_4,omnes,matched_tilde,phases, 2, epsilon, validity, cutoff, SubtractionConstant::a1, max_subs),
//C-violating
disp_g1_1(mass_1,mass_4,omnes,matched_tilde,phases, 1, epsilon, validity, cutoff, SubtractionConstant::g1, max_subs),
disp_h1_1(mass_1,mass_4,omnes,matched_tilde,phases, 1, epsilon, validity, cutoff, SubtractionConstant::h1, max_subs),
disp_h1_2(mass_1,mass_4,omnes,matched_tilde,phases, 2, epsilon, validity, cutoff, SubtractionConstant::h1, max_subs)
{}



BasisAmplitudeEtapEtaPiPi::BasisAmplitudeEtapEtaPiPi(double mass_1, double mass_3, double mass_4,
                                FunctionSetOmnes omnes,
                                FunctionSet matched_tilde,
                                phase::PhaseEtapEtaPiPi phases,
                                double cutoff,double epsilon, double validity, int max_subs)
:
BasisAmplitude(mass_1,mass_1,mass_3,mass_4,omnes,matched_tilde,cutoff,epsilon,validity,max_subs),
gamma(mass_1,mass_1,mass_3,mass_4),
/// Notation for the last three digits: e.g. ang_a0_XYZ has total isospin X, isospin of two-particle intermediate state Y and partial wave Z
disp_a0_000(mass_1,mass_3,mass_4,omnes,matched_tilde,phases, 000, epsilon,validity,cutoff, SubtractionConstant::sub_a0, max_subs),
disp_a0_010(mass_1,mass_3,mass_4,omnes,matched_tilde,phases, 010, epsilon,validity,cutoff, SubtractionConstant::sub_a0, max_subs),
disp_b0_000(mass_1,mass_3,mass_4,omnes,matched_tilde,phases, 000, epsilon,validity,cutoff, SubtractionConstant::sub_b0, max_subs),
disp_b0_010(mass_1,mass_3,mass_4,omnes,matched_tilde,phases, 010, epsilon,validity,cutoff, SubtractionConstant::sub_b0, max_subs),
disp_c0_000(mass_1,mass_3,mass_4,omnes,matched_tilde,phases, 000, epsilon,validity,cutoff, SubtractionConstant::sub_c0, max_subs),
disp_c0_010(mass_1,mass_3,mass_4,omnes,matched_tilde,phases, 010, epsilon,validity,cutoff, SubtractionConstant::sub_c0, max_subs),
disp_d0_000(mass_1,mass_3,mass_4,omnes,matched_tilde,phases, 000, epsilon,validity,cutoff, SubtractionConstant::sub_d0, max_subs),
disp_d0_010(mass_1,mass_3,mass_4,omnes,matched_tilde,phases, 010, epsilon,validity,cutoff, SubtractionConstant::sub_d0, max_subs),
// C-violating
disp_a1_111(mass_1,mass_3,mass_4,omnes,matched_tilde,phases, 111, epsilon,validity,cutoff, SubtractionConstant::sub_a1, max_subs),
disp_a1_110(mass_1,mass_3,mass_4,omnes,matched_tilde,phases, 110, epsilon,validity,cutoff, SubtractionConstant::sub_a1, max_subs),
disp_b1_111(mass_1,mass_3,mass_4,omnes,matched_tilde,phases, 111, epsilon,validity,cutoff, SubtractionConstant::sub_b1, max_subs),
disp_b1_110(mass_1,mass_3,mass_4,omnes,matched_tilde,phases, 110, epsilon,validity,cutoff, SubtractionConstant::sub_b1, max_subs)
{}

BasisAmplitudeV3Pi::BasisAmplitudeV3Pi(double mass_1, double mass_4,
                               FunctionSetOmnes omnes,
                               FunctionSet matched_tilde,
                               phase::PhaseV3Pi phases,
                               double cutoff,double epsilon, double validity, int max_subs)
:
BasisAmplitude(mass_1,mass_1,mass_1,mass_4,omnes,matched_tilde,cutoff,epsilon,validity,max_subs),
polar_egg(mass_4, mass_1),
disp_a1(mass_1,mass_4,omnes,matched_tilde,phases, 1, epsilon,validity,cutoff, SubtractionConstant::a1, max_subs),
disp_b1(mass_1,mass_4,omnes,matched_tilde,phases, 1, epsilon,validity,cutoff, SubtractionConstant::b1, max_subs)
{}

BasisAmplitudeX3Pi::BasisAmplitudeX3Pi(double mass_1, double mass_4,
                               FunctionSetOmnes omnes,
                               FunctionSet matched_tilde,
                               phase::PhaseX3Pi phases,
                               double cutoff,double epsilon, double validity, int max_subs)
:
BasisAmplitude(mass_1,mass_1,mass_1,mass_4,omnes,matched_tilde,cutoff,epsilon,validity,max_subs),
polar_egg(mass_4, mass_1),
disp_a1(mass_1,mass_4,omnes,matched_tilde,phases, 1, epsilon,validity,cutoff, SubtractionConstant::a1, max_subs),
disp_b1(mass_1,mass_4,omnes,matched_tilde,phases, 1, epsilon,validity,cutoff, SubtractionConstant::b1, max_subs)
{}

BasisAmplitudeT3Pi::BasisAmplitudeT3Pi(double mass_1, double mass_4,
                               FunctionSetOmnes omnes,
                               FunctionSet matched_tilde,
                               phase::PhaseT3Pi phases,
                               double cutoff,double epsilon, double validity, int max_subs)
:
BasisAmplitude(mass_1,mass_1,mass_1,mass_4,omnes,matched_tilde,cutoff,epsilon,validity,max_subs),
polar_egg(mass_4, mass_1),
disp_a1(mass_1,mass_4,omnes,matched_tilde,phases, 1, epsilon,validity,cutoff, SubtractionConstant::a1, max_subs),
disp_b1(mass_1,mass_4,omnes,matched_tilde,phases, 1, epsilon,validity,cutoff, SubtractionConstant::b1, max_subs)
{}


//-------------------------------------------------------------------------------------------
//-- Homogeneous Solution--------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-- ADJUST FOR DIFFERENT SUBTRACTION SCHEME ------------------------------------------------

Complex HomogeneousSolutionEta3Pi::initial_guess(double s, int isospin, SubtractionConstant sub,
                                           Setting evaluation)const{
    switch (sub) {
        // Standard Model (total isospin 1, C-conserving)
        case SubtractionConstant::a0:
            if (isospin==0) {return omnes(s, isospin, evaluation);}
            if (isospin==1) {return 0.;}
            if (isospin==2) {return 0.;}
            throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};
            break;
            
        case SubtractionConstant::b0:
            if (isospin==0) {
                if (evaluation==Setting::egg) {return omnes(s, isospin, evaluation)* polar_egg(s);}
                return s*omnes(s, isospin, evaluation);
            }
            if (isospin==1) {return 0.;}
            if (isospin==2) {return 0.;}
            throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};
            break;
            
        case SubtractionConstant::a1:
            if (isospin==0) {return 0.;}
            if (isospin==1) {return omnes(s, isospin, evaluation);}
            if (isospin==2) {return 0.;}
            throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};
            break;
        
        // Beyond Standard Model (total isospin 0, C-violating)
        case SubtractionConstant::g1:
            if (isospin==1) {return omnes(s, isospin, evaluation);}
            throw std::domain_error{"For total isospin 0 the C-violating amplitude has to have isospin 1."};
            break;
            
        // Beyond Standard Model (total isospin 2, C-violating)
        case SubtractionConstant::h1:
            if (isospin==1) {return omnes(s, isospin, evaluation);}
            if (isospin==2) {return 0.;}
            throw std::domain_error{"For total isospin 2 the C-violating amplitude has to have isospin 1 or 2."};
            break;
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
            break;
    }
}

Complex HomogeneousSolutionEtap3Pi::initial_guess(double s, int isospin, SubtractionConstant sub,
                                           Setting evaluation)const{
    switch (sub) {
        // Standard Model (total isospin 1, C-conserving)
        case SubtractionConstant::a0:
            if (isospin==0) {return omnes(s, isospin, evaluation);}
            if (isospin==1) {return 0.;}
            if (isospin==2) {return 0.;}
            throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};
            break;
            
        case SubtractionConstant::b0:
            if (isospin==0) {
                if (evaluation==Setting::egg) {return omnes(s, isospin, evaluation)* polar_egg(s);}
                return s*omnes(s, isospin, evaluation);
            }
            if (isospin==1) {return 0.;}
            if (isospin==2) {return 0.;}
            throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};
            break;
            
        case SubtractionConstant::a1:
            if (isospin==0) {return 0.;}
            if (isospin==1) {return omnes(s, isospin, evaluation);}
            if (isospin==2) {return 0.;}
            throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};
            break;
        
        // Beyond Standard Model (total isospin 0, C-violating)
        case SubtractionConstant::g1:
            if (isospin==1) {return omnes(s, isospin, evaluation);}
            throw std::domain_error{"For total isospin 0 the C-violating amplitude has to have isospin 1."};
            break;
            
        // Beyond Standard Model (total isospin 2, C-violating)
        case SubtractionConstant::h1:
            if (isospin==1) {return omnes(s, isospin, evaluation);}
            if (isospin==2) {return 0.;}
            throw std::domain_error{"For total isospin 2 the C-violating amplitude has to have isospin 1 or 2."};
            break;
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
            break;
    }
}




Complex HomogeneousSolutionEtapEtaPiPi::initial_guess(double s, int isospin, SubtractionConstant sub,
                                           Setting evaluation)const{
    switch (sub) {
        // Standard Model (total isospin 1, C-conserving)
        case SubtractionConstant::sub_a0:
            if (isospin==000) {return omnes(s, isospin, evaluation);}
            if (isospin==010) {return 0.;}
            throw std::domain_error{"Wrong isospin for homogeneous solution."};
            break;
            
        case SubtractionConstant::sub_b0:
            if (isospin==000) {
                switch (evaluation) {
                    case Setting::above: case Setting::below: return s * omnes(s, isospin, evaluation); break;
                    case Setting::egg:                        throw  std::domain_error{"Demanded evaluation is not needed for this process."}; break;
                    default:                                  return gamma(s, evaluation) * omnes(s, isospin, evaluation); break;
                }
            }
            if (isospin==010) {return 0.;}
            throw std::domain_error{"Wrong isospin for homogeneous solution."};
            break;
            
        case SubtractionConstant::sub_c0:
            if (isospin==000) {
                switch (evaluation) {
                    case Setting::above: case Setting::below: return std::pow(s,2.) * omnes(s, isospin, evaluation); break;
                    case Setting::egg:                        throw  std::domain_error{"Demanded evaluation is not needed for this process."}; break;
                    default:                                  return std::pow(gamma(s, evaluation),2.) * omnes(s, isospin, evaluation); break;
                }
            }
            if (isospin==010) {return 0.;}
            throw std::domain_error{"Wrong isospin for homogeneous solution."};
            break;
            
            
        case SubtractionConstant::sub_d0:
            if (isospin==000) {return 0.;}
            if (isospin==010) {
                switch (evaluation) {
                    case Setting::above: case Setting::below: return std::pow(s,2.) * omnes(s, isospin, evaluation); break;
                    case Setting::egg:                        throw  std::domain_error{"Demanded evaluation is not needed for this process."}; break;
                    default:                                  return std::pow(gamma(s, evaluation),2.) * omnes(s, isospin, evaluation); break;
                }
            }
            throw std::domain_error{"Wrong isospin for homogeneous solution."};
            break;
            
            
            // Beyond Standard Model (total isospin 1)
        case SubtractionConstant::sub_a1:
            if (isospin==111) {return omnes(s, isospin, evaluation);}
            if (isospin==110) {return 0.;}
            throw std::domain_error{"Wrong isospin for homogeneous solution."};
            break;
            
            
        case SubtractionConstant::sub_b1:
            if (isospin==111) {return 0.;}
            if (isospin==110) {
                switch (evaluation) {
                    case Setting::above: case Setting::below: return std::pow(s,1.) * omnes(s, isospin, evaluation); break;
                    case Setting::egg:                        throw  std::domain_error{"Demanded evaluation is not needed for this process."}; break;
                    default:                                  return std::pow(gamma(s, evaluation),1.) * omnes(s, isospin, evaluation); break;
                }
            }
            throw std::domain_error{"Wrong isospin for homogeneous solution."};
            break;
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
            break;
    }// switch
}

Complex HomogeneousSolutionV3Pi::initial_guess(double s, int isospin, SubtractionConstant sub,
                                           Setting evaluation)const{
    switch (sub) {
        // Standard Model (total isospin 1, C-conserving)
        case SubtractionConstant::a1:
            return omnes(s, isospin, evaluation);
            
        case SubtractionConstant::b1:
            if (evaluation==Setting::egg) {return omnes(s, isospin, evaluation)* polar_egg(s);}
            return s*omnes(s, isospin, evaluation);
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
            break;
    }
}

Complex HomogeneousSolutionX3Pi::initial_guess(double s, int isospin, SubtractionConstant sub,
                                           Setting evaluation)const{
    switch (sub) {
        // Standard Model (total isospin 1, C-conserving)
        case SubtractionConstant::a1:
            return omnes(s, isospin, evaluation);
            
        case SubtractionConstant::b1:
            if (evaluation==Setting::egg) {return omnes(s, isospin, evaluation)* polar_egg(s);}
            return s*omnes(s, isospin, evaluation);
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
            break;
    }
}

Complex HomogeneousSolutionT3Pi::initial_guess(double s, int isospin, SubtractionConstant sub,
                                           Setting evaluation)const{
    switch (sub) {
        // Standard Model (total isospin 1, C-conserving)
        case SubtractionConstant::a1:
            return omnes(s, isospin, evaluation);
            
        case SubtractionConstant::b1:
            if (evaluation==Setting::egg) {return omnes(s, isospin, evaluation)* polar_egg(s);}
            return s*omnes(s, isospin, evaluation);
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
            break;
    }
}
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------





//-------------------------------------------------------------------------------------------
//-- Tilde Function -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------


Complex BasisTildeEta3Pi::matched_tilde(double s, int isospin, SubtractionConstant sub)const{
    // Standard Model (total isospin 1, C-conserving)
    switch (sub) {
        case SubtractionConstant::a0:
            if (isospin==0) {return ang_a0_0(s);}
            if (isospin==1) {return ang_a0_1(s);}
            if (isospin==2) {return ang_a0_2(s);}
            else{throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};};
            
        case SubtractionConstant::b0:
            if (isospin==0) {return ang_b0_0(s);}
            if (isospin==1) {return ang_b0_1(s);}
            if (isospin==2) {return ang_b0_2(s);}
            else{throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};};
            
        case SubtractionConstant::a1:
            if (isospin==0) {return ang_a1_0(s);}
            if (isospin==1) {return ang_a1_1(s);}
            if (isospin==2) {return ang_a1_2(s);}
            else{throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};};
            
        // Beyond Standard Model (total isospin 0, C-violating)
        case SubtractionConstant::g1:
            if (isospin==1) {return ang_g1_1(s);}
            else{throw std::domain_error{"For total isospin 0 the C-violating amplitude has to have isospin 1."};};
            
        // Beyond Standard Model (total isospin 2, C-violating)
        case SubtractionConstant::h1:
            if (isospin==1) {return ang_h1_1(s);}
            if (isospin==2) {return ang_h1_2(s);}
            else{throw std::domain_error{"For total isospin 2 the C-violating amplitude has to have isospin 1 or 2."};};
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
    }//switch
}

Complex BasisTildeEtap3Pi::matched_tilde(double s, int isospin, SubtractionConstant sub)const{
    // Standard Model (total isospin 1, C-conserving)
    switch (sub) {
        case SubtractionConstant::a0:
            if (isospin==0) {return ang_a0_0(s);}
            if (isospin==1) {return ang_a0_1(s);}
            if (isospin==2) {return ang_a0_2(s);}
            else{throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};};
            
        case SubtractionConstant::b0:
            if (isospin==0) {return ang_b0_0(s);}
            if (isospin==1) {return ang_b0_1(s);}
            if (isospin==2) {return ang_b0_2(s);}
            else{throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};};
            
        case SubtractionConstant::a1:
            if (isospin==0) {return ang_a1_0(s);}
            if (isospin==1) {return ang_a1_1(s);}
            if (isospin==2) {return ang_a1_2(s);}
            else{throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};};
            
        // Beyond Standard Model (total isospin 0, C-violating)
        case SubtractionConstant::g1:
            if (isospin==1) {return ang_g1_1(s);}
            else{throw std::domain_error{"For total isospin 0 the C-violating amplitude has to have isospin 1."};};
            
        // Beyond Standard Model (total isospin 2, C-violating)
        case SubtractionConstant::h1:
            if (isospin==1) {return ang_h1_1(s);}
            if (isospin==2) {return ang_h1_2(s);}
            else{throw std::domain_error{"For total isospin 2 the C-violating amplitude has to have isospin 1 or 2."};};
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
    }//switch
}


Complex BasisTildeEtapEtaPiPi::matched_tilde(double s, int isospin, SubtractionConstant sub)const{    
    // Standard Model (total isospin 1, C-conserving)
    switch (sub) {
        case SubtractionConstant::sub_a0:
            if (isospin==000) {return ang_a0_000(s);}
            if (isospin==010) {return ang_a0_010(s);}
            else{throw std::domain_error{"Wrong isospin."};};
            
        case SubtractionConstant::sub_b0:
            if (isospin==000) {return ang_b0_000(s);}
            if (isospin==010) {return ang_b0_010(s);}
            else{throw std::domain_error{"Wrong isospin."};};
            
        case SubtractionConstant::sub_c0:
            if (isospin==000) {return ang_c0_000(s);}
            if (isospin==010) {return ang_c0_010(s);}
            else{throw std::domain_error{"Wrong isospin."};};
            
        case SubtractionConstant::sub_d0:
            if (isospin==000) {return ang_d0_000(s);}
            if (isospin==010) {return ang_d0_010(s);}
            else{throw std::domain_error{"Wrong isospin."};};
            
        // Beyond Standard Model (total isospin 1, C-violating)
        case SubtractionConstant::sub_a1:
            if (isospin==111) {return ang_a1_111(s);}
            if (isospin==110) {return ang_a1_110(s);}
            else{throw std::domain_error{"Wrong isospin."};};
            
        case SubtractionConstant::sub_b1:
            if (isospin==111) {return ang_b1_111(s);}
            if (isospin==110) {return ang_b1_110(s);}
            else{throw std::domain_error{"Wrong isospin."};};
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
    }//switch
}

Complex BasisTildeV3Pi::matched_tilde(double s, int isospin, SubtractionConstant sub)const{
    // Standard Model (total isospin 1, C-conserving)
    switch (sub) {
        case SubtractionConstant::a1:
            return ang_a1(s);
            
        case SubtractionConstant::b1:
            return ang_b1(s);

        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
    }//switch
}

Complex BasisTildeX3Pi::matched_tilde(double s, int isospin, SubtractionConstant sub)const{
    // Standard Model (total isospin 1, C-conserving)
    switch (sub) {
        case SubtractionConstant::a1:
            return ang_a1(s);
            
        case SubtractionConstant::b1:
            return ang_b1(s);

        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
    }//switch
}

Complex BasisTildeT3Pi::matched_tilde(double s, int isospin, SubtractionConstant sub)const{
    // Standard Model (total isospin 1, C-conserving)
    switch (sub) {
        case SubtractionConstant::a1:
            return ang_a1(s);
            
        case SubtractionConstant::b1:
            return ang_b1(s);

        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
    }//switch
}
//-------------------------------------------------------------------------
//-------------------------------------------------------------------------





//-------------------------------------------------------------------------------------------
//-- BasisAmplitude -------------------------------------------------------------------------
//-------------------------------------------------------------------------------------------
//-- ADJUST FOR DIFFERENT SUBTRACTION SCHEME ------------------------------------------------

// To evaluate along the complex egg-like contour, one has to express the basis amplitude for this
// case in terms of the real-valued curve-parameter.
// Note that this curve parameter is not 's' anymore but rather an angle in radian.

Complex BasisAmplitudeEta3Pi::dispersion_integral(double s, int isospin, SubtractionConstant sub,
                                            Setting evaluation)const{
    switch (sub) {
        // Standard Model (total isospin 1, C-conserving)
        case SubtractionConstant::a0:
            if (isospin==0) {return disp_a0_0(s, evaluation);}
            if (isospin==1) {return disp_a0_1(s, evaluation);}
            if (isospin==2) {return disp_a0_2(s, evaluation);}
            else{throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};};
            
        case SubtractionConstant::b0:
            if (isospin==0) {return disp_b0_0(s, evaluation);}
            if (isospin==1) {return disp_b0_1(s, evaluation);}
            if (isospin==2) {return disp_b0_2(s, evaluation);}
            else{throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};};
            
        case SubtractionConstant::a1:
            if (isospin==0) {return disp_a1_0(s, evaluation);}
            if (isospin==1) {return disp_a1_1(s, evaluation);}
            if (isospin==2) {return disp_a1_2(s, evaluation);}
            else{throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};};
            
        // Beyond Standard Model (total isospin 0, C-violating)
        case SubtractionConstant::g1:
            if (isospin==1) {return disp_g1_1(s, evaluation);}
            else{throw std::domain_error{"For total isospin 0 the C-violating amplitude has to have isospin 1."};};
            
        // Beyond Standard Model (total isospin 2, C-violating)
        case SubtractionConstant::h1:
            if (isospin==1) {return disp_h1_1(s, evaluation);}
            if (isospin==2) {return disp_h1_2(s, evaluation);}
            else{throw std::domain_error{"For total isospin 2 the C-violating amplitude has to have isospin 1 or 2."};};
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
    }//switch
}

Complex BasisAmplitudeEtap3Pi::dispersion_integral(double s, int isospin, SubtractionConstant sub,
                                            Setting evaluation)const{
    switch (sub) {
        // Standard Model (total isospin 1, C-conserving)
        case SubtractionConstant::a0:
            if (isospin==0) {return disp_a0_0(s, evaluation);}
            if (isospin==1) {return disp_a0_1(s, evaluation);}
            if (isospin==2) {return disp_a0_2(s, evaluation);}
            else{throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};};
            
        case SubtractionConstant::b0:
            if (isospin==0) {return disp_b0_0(s, evaluation);}
            if (isospin==1) {return disp_b0_1(s, evaluation);}
            if (isospin==2) {return disp_b0_2(s, evaluation);}
            else{throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};};
            
        case SubtractionConstant::a1:
            if (isospin==0) {return disp_a1_0(s, evaluation);}
            if (isospin==1) {return disp_a1_1(s, evaluation);}
            if (isospin==2) {return disp_a1_2(s, evaluation);}
            else{throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};};
            
        // Beyond Standard Model (total isospin 0, C-violating)
        case SubtractionConstant::g1:
            if (isospin==1) {return disp_g1_1(s, evaluation);}
            else{throw std::domain_error{"For total isospin 0 the C-violating amplitude has to have isospin 1."};};
            
        // Beyond Standard Model (total isospin 2, C-violating)
        case SubtractionConstant::h1:
            if (isospin==1) {return disp_h1_1(s, evaluation);}
            if (isospin==2) {return disp_h1_2(s, evaluation);}
            else{throw std::domain_error{"For total isospin 2 the C-violating amplitude has to have isospin 1 or 2."};};
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
    }//switch
}


Complex BasisAmplitudeEtapEtaPiPi::dispersion_integral(double s, int isospin, SubtractionConstant sub,
                                            Setting evaluation)const{
    switch (sub) {
        // Standard Model (total isospin 1, C-conserving)
        case SubtractionConstant::sub_a0:
            if (isospin==000) {return disp_a0_000(s, evaluation);}
            if (isospin==010) {return disp_a0_010(s, evaluation);}
            else{throw std::domain_error{"Wrong isospin."};};
            
        case SubtractionConstant::sub_b0:
            if (isospin==000) {return disp_b0_000(s, evaluation);}
            if (isospin==010) {return disp_b0_010(s, evaluation);}
            else{throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};};
            
        case SubtractionConstant::sub_c0:
            if (isospin==000) {return disp_c0_000(s, evaluation);}
            if (isospin==010) {return disp_c0_010(s, evaluation);}
            else{throw std::domain_error{"Wrong isospin."};};
            
        case SubtractionConstant::sub_d0:
            if (isospin==000) {return disp_d0_000(s, evaluation);}
            if (isospin==010) {return disp_d0_010(s, evaluation);}
            else{throw std::domain_error{"Wrong isospin."};};
            
        // Beyond Standard Model (total isospin 1, C-violating)
        case SubtractionConstant::sub_a1:
            if (isospin==111) {return disp_a1_111(s, evaluation);}
            if (isospin==110) {return disp_a1_110(s, evaluation);}
            else{throw std::domain_error{"Wrong isospin."};};
            
        case SubtractionConstant::sub_b1:
            if (isospin==111) {return disp_b1_111(s, evaluation);}
            if (isospin==110) {return disp_b1_110(s, evaluation);}
            else{throw std::domain_error{"Wrong isospin."};};
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
    }//switch
}

Complex BasisAmplitudeV3Pi::dispersion_integral(double s, int isospin, SubtractionConstant sub,
                                            Setting evaluation)const{
    switch (sub) {
        // Standard Model (total isospin 1, C-conserving)
        case SubtractionConstant::a1:
            return disp_a1(s, evaluation);

        case SubtractionConstant::b1:
            return disp_b1(s, evaluation);
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
    }//switch
}

Complex BasisAmplitudeX3Pi::dispersion_integral(double s, int isospin, SubtractionConstant sub,
                                            Setting evaluation)const{
    switch (sub) {
        // Standard Model (total isospin 1, C-conserving)
        case SubtractionConstant::a1:
            return disp_a1(s, evaluation);

        case SubtractionConstant::b1:
            return disp_b1(s, evaluation);
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
    }//switch
}

Complex BasisAmplitudeT3Pi::dispersion_integral(double s, int isospin, SubtractionConstant sub,
                                            Setting evaluation)const{
    switch (sub) {
        // Standard Model (total isospin 1, C-conserving)
        case SubtractionConstant::a1:
            return disp_a1(s, evaluation);

        case SubtractionConstant::b1:
            return disp_b1(s, evaluation);
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
    }//switch
}



Complex BasisAmplitudeEta3Pi::helper(double s, int isospin, SubtractionConstant sub, Setting evaluation,
                                              double number_subtractions, double power_of_s,
                                              double polynomial_on_off)const{
    
    if (polynomial_on_off!=1 && polynomial_on_off!=0) {
            throw std::domain_error{"Input for 'polynomial_on_off' not allowed."};
    }
    
    switch (evaluation) {
        case Setting::above: case Setting::below:
            return
            omnes(s,isospin,evaluation)* (std::pow(s,power_of_s) *polynomial_on_off
                                          + std::pow(s,number_subtractions)/pi()
                                          * dispersion_integral(s,isospin,sub,evaluation));
            
        case Setting::egg:
            return
            omnes(s,isospin,evaluation)* (std::pow(polar_egg(s),power_of_s) *polynomial_on_off
                                          + std::pow(polar_egg(s),number_subtractions)/pi()
                                          * dispersion_integral(s,isospin,sub,evaluation));
            
        default:
            throw std::domain_error{"Chosen evaluation not needed for the desired process."};
    }
}

Complex BasisAmplitudeEtap3Pi::helper(double s, int isospin, SubtractionConstant sub, Setting evaluation,
                                              double number_subtractions, double power_of_s,
                                              double polynomial_on_off)const{
    
    if (polynomial_on_off!=1 && polynomial_on_off!=0) {
            throw std::domain_error{"Input for 'polynomial_on_off' not allowed."};
    }
    
    switch (evaluation) {
        case Setting::above: case Setting::below:
            return
            omnes(s,isospin,evaluation)* (std::pow(s,power_of_s) *polynomial_on_off
                                          + std::pow(s,number_subtractions)/pi()
                                          * dispersion_integral(s,isospin,sub,evaluation));
            
        case Setting::egg:
            return
            omnes(s,isospin,evaluation)* (std::pow(polar_egg(s),power_of_s) *polynomial_on_off
                                          + std::pow(polar_egg(s),number_subtractions)/pi()
                                          * dispersion_integral(s,isospin,sub,evaluation));
            
        default:
            throw std::domain_error{"Chosen evaluation not needed for the desired process."};
    }
}

Complex BasisAmplitudeEtapEtaPiPi::helper(double s, int isospin, SubtractionConstant sub, Setting evaluation,
                                              double number_subtractions, double power_of_s,
                                              double polynomial_on_off)const{
    
    if (polynomial_on_off!=1 && polynomial_on_off!=0) {
            throw std::domain_error{"Input for 'polynomial_on_off' not allowed."};
    }
    
    switch (evaluation) {
        case Setting::above: case Setting::below:
            return
            omnes(s,isospin,evaluation)* (std::pow(s,power_of_s) *polynomial_on_off
                                          + std::pow(s,number_subtractions)/pi()
                                          * dispersion_integral(s,isospin,sub,evaluation));
            
        case Setting::egg:
            throw std::domain_error{"Chosen evaluation not needed for the desired process."};
            
        default:
            return
            omnes(s,isospin,evaluation)* (std::pow(gamma(s, evaluation),power_of_s) *polynomial_on_off
                                          + std::pow(gamma(s, evaluation),number_subtractions)/pi()
                                          * dispersion_integral(s,isospin,sub,evaluation));
    }
}

Complex BasisAmplitudeV3Pi::helper(double s, int isospin, SubtractionConstant sub, Setting evaluation,
                                              double number_subtractions, double power_of_s,
                                              double polynomial_on_off)const{
    
    if (polynomial_on_off!=1 && polynomial_on_off!=0) {
            throw std::domain_error{"Input for 'polynomial_on_off' not allowed."};
    }
    
    switch (evaluation) {
        case Setting::above: case Setting::below:
            return
            omnes(s,isospin,evaluation)* (std::pow(s,power_of_s) *polynomial_on_off
                                          + std::pow(s,number_subtractions)/pi()
                                          * dispersion_integral(s,isospin,sub,evaluation));
            
        case Setting::egg:
            return
            omnes(s,isospin,evaluation)* (std::pow(polar_egg(s),power_of_s) *polynomial_on_off
                                          + std::pow(polar_egg(s),number_subtractions)/pi()
                                          * dispersion_integral(s,isospin,sub,evaluation));
            
        default:
            throw std::domain_error{"Chosen evaluation not needed for the desired process."};
    }
}

Complex BasisAmplitudeX3Pi::helper(double s, int isospin, SubtractionConstant sub, Setting evaluation,
                                              double number_subtractions, double power_of_s,
                                              double polynomial_on_off)const{
    
    if (polynomial_on_off!=1 && polynomial_on_off!=0) {
            throw std::domain_error{"Input for 'polynomial_on_off' not allowed."};
    }
    
    switch (evaluation) {
        case Setting::above: case Setting::below:
            return
            omnes(s,isospin,evaluation)* (std::pow(s,power_of_s) *polynomial_on_off
                                          + std::pow(s,number_subtractions)/pi()
                                          * dispersion_integral(s,isospin,sub,evaluation));
            
        case Setting::egg:
            return
            omnes(s,isospin,evaluation)* (std::pow(polar_egg(s),power_of_s) *polynomial_on_off
                                          + std::pow(polar_egg(s),number_subtractions)/pi()
                                          * dispersion_integral(s,isospin,sub,evaluation));
            
        default:
            throw std::domain_error{"Chosen evaluation not needed for the desired process."};
    }
}

Complex BasisAmplitudeT3Pi::helper(double s, int isospin, SubtractionConstant sub, Setting evaluation,
                                              double number_subtractions, double power_of_s,
                                              double polynomial_on_off)const{
    
    if (polynomial_on_off!=1 && polynomial_on_off!=0) {
            throw std::domain_error{"Input for 'polynomial_on_off' not allowed."};
    }
    
    switch (evaluation) {
        case Setting::above: case Setting::below:
            return
            omnes(s,isospin,evaluation)* (std::pow(s,power_of_s) *polynomial_on_off
                                          + std::pow(s,number_subtractions)/pi()
                                          * dispersion_integral(s,isospin,sub,evaluation));
            
        case Setting::egg:
            return
            omnes(s,isospin,evaluation)* (std::pow(polar_egg(s),power_of_s) *polynomial_on_off
                                          + std::pow(polar_egg(s),number_subtractions)/pi()
                                          * dispersion_integral(s,isospin,sub,evaluation));
            
        default:
            throw std::domain_error{"Chosen evaluation not needed for the desired process."};
    }
}





Complex BasisAmplitudeEta3Pi::amplitude(double s, int isospin, SubtractionConstant sub,
                                  Setting evaluation)const{
    
    switch (sub) {
        // Standard Model (total isospin 1, C-conserving)
        case SubtractionConstant::a0:
            switch (isospin) {
                case 0: return helper(s, isospin, sub, evaluation, 2., 0., 1); 
                case 1: case 2: return helper(s, isospin, sub, evaluation, 1.); 
                default:
                    throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};
            }
            
        case SubtractionConstant::b0:
            switch (isospin) {
                case 0: return helper(s, isospin, sub, evaluation, 2., 1., 1); 
                case 1: case 2: return helper(s, isospin, sub, evaluation, 1.); 
                default:
                    throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};
            }
            
        case SubtractionConstant::a1:
            switch (isospin) {
                case 0: return helper(s, isospin, sub, evaluation, 2.);
                case 1: return helper(s, isospin, sub, evaluation, 1., 0., 1); 
                case 2: return helper(s, isospin, sub, evaluation, 1.); 
                default:
                    throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};
            }
        
        // Beyond Standard Model (total isospin 0, C-violating)
        case SubtractionConstant::g1:
            switch (isospin) {
                case 1: return helper(s, isospin, sub, evaluation, 1., 0., 1); 
                default:
                    throw std::domain_error{"For total isospin 0 the C-violating amplitude has to have isospin 1."};
            }
            
        // Beyond Standard Model (total isospin 2, C-violating)
        case SubtractionConstant::h1:
            switch (isospin) {
                case 1: return helper(s, isospin, sub, evaluation, 1., 0., 1); 
                case 2: return helper(s, isospin, sub, evaluation, 1.); 
                default:
                    throw std::domain_error{"For total isospin 2 the C-violating amplitude has to have isospin 1 or 2."};
            }
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
    }//switch
}

Complex BasisAmplitudeEtap3Pi::amplitude(double s, int isospin, SubtractionConstant sub,
                                  Setting evaluation)const{
    
    switch (sub) {
        // Standard Model (total isospin 1, C-conserving)
        case SubtractionConstant::a0:
            switch (isospin) {
                case 0: return helper(s, isospin, sub, evaluation, 2., 0., 1); 
                case 1: case 2: return helper(s, isospin, sub, evaluation, 1.); 
                default:
                    throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};
            }
            
        case SubtractionConstant::b0:
            switch (isospin) {
                case 0: return helper(s, isospin, sub, evaluation, 2., 1., 1); 
                case 1: case 2: return helper(s, isospin, sub, evaluation, 1.); 
                default:
                    throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};

            }
            
        case SubtractionConstant::a1:
            switch (isospin) {
                case 0: return helper(s, isospin, sub, evaluation, 2.); 
                case 1: return helper(s, isospin, sub, evaluation, 1., 0., 1); 
                case 2: return helper(s, isospin, sub, evaluation, 1.); 
                default:
                    throw std::domain_error{"The two-pion state must have either isospin 0, 1 or 2"};
            }
        
        // Beyond Standard Model (total isospin 0, C-violating)
        case SubtractionConstant::g1:
            switch (isospin) {
                case 1: return helper(s, isospin, sub, evaluation, 1., 0., 1); 
                default:
                    throw std::domain_error{"For total isospin 0 the C-violating amplitude has to have isospin 1."};
            }
            
        // Beyond Standard Model (total isospin 2, C-violating)
        case SubtractionConstant::h1:
            switch (isospin) {
                case 1: return helper(s, isospin, sub, evaluation, 1., 0., 1); 
                case 2: return helper(s, isospin, sub, evaluation, 1.); 
                default:
                    throw std::domain_error{"For total isospin 2 the C-violating amplitude has to have isospin 1 or 2."};
            }
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
    }//switch
}

Complex BasisAmplitudeEtapEtaPiPi::amplitude(double s, int isospin, SubtractionConstant sub,
                                  Setting evaluation)const{
    
    switch (sub) {
        // Standard Model (total isospin 1, C-conserving)
        case SubtractionConstant::sub_a0:
            switch (isospin) {
                case 000:
                    return helper(s,isospin,sub,evaluation, 3., 0., 1);
                    
                case 010:
                    return helper(s,isospin,sub,evaluation, 3.);
            
                default:
                    throw std::domain_error{"Isospin does not contribute to the chosen basis function."};
            }
            
        case SubtractionConstant::sub_b0:
            switch (isospin) {
                case 000:
                    return helper(s,isospin,sub,evaluation, 3., 1., 1);
                    
                case 010:
                    return helper(s,isospin,sub,evaluation, 3.);
            
                default:
                    throw std::domain_error{"Isospin does not contribute to the chosen basis function."};
            }
            
        case SubtractionConstant::sub_c0:
            switch (isospin) {
                case 000:
                    return helper(s,isospin,sub,evaluation, 3., 2., 1);
                    
                case 010:
                    return helper(s,isospin,sub,evaluation, 3.);
            
                default:
                    throw std::domain_error{"Isospin does not contribute to the chosen basis function."};
            }
            
        case SubtractionConstant::sub_d0:
            switch (isospin) {
                case 000:
                    return helper(s,isospin,sub,evaluation, 3.);
                    
                case 010:
                    return helper(s,isospin,sub,evaluation, 3., 2., 1);
            
                default:
                    throw std::domain_error{"Isospin does not contribute to the chosen basis function."};
            }
            
            
        case SubtractionConstant::sub_a1:
            switch (isospin) {
                case 111:
                    return helper(s,isospin,sub,evaluation, 1., 0., 1);
                    
                case 110:
                    return helper(s,isospin,sub,evaluation, 2.);
                              
                default:
                    throw std::domain_error{"Isospin does not contribute to the chosen basis function."};                  
            }
            
        case SubtractionConstant::sub_b1:
            switch (isospin) {
                case 111:
                    return helper(s,isospin,sub,evaluation, 1.);
                                       
                case 110:
                    return helper(s,isospin,sub,evaluation, 2., 1., 1);                
            
                default:
                    throw std::domain_error{"Isospin does not contribute to the chosen basis function."};                 
            }
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
            
    }//switch
}

Complex BasisAmplitudeV3Pi::amplitude(double s, int isospin, SubtractionConstant sub,
                                  Setting evaluation)const{
    
    switch (sub) {
        // Standard Model (total isospin 1)
        case SubtractionConstant::a1:
            return helper(s, isospin, sub, evaluation, max_subs, 0., 1);
            
        case SubtractionConstant::b1:
            return helper(s, isospin, sub, evaluation, max_subs, 1., 1);
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
            
    }//switch
}

Complex BasisAmplitudeX3Pi::amplitude(double s, int isospin, SubtractionConstant sub,
                                  Setting evaluation)const{
    
    switch (sub) {
        // Standard Model (total isospin 1)
        case SubtractionConstant::a1:
            return helper(s, isospin, sub, evaluation, max_subs, 0., 1);
            
        case SubtractionConstant::b1:
            return helper(s, isospin, sub, evaluation, max_subs, 1., 1);
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
            
    }//switch
}

Complex BasisAmplitudeT3Pi::amplitude(double s, int isospin, SubtractionConstant sub,
                                  Setting evaluation)const{
    
    switch (sub) {
        // Standard Model (total isospin 1)
        case SubtractionConstant::a1:
            return helper(s, isospin, sub, evaluation, max_subs, 0., 1);
            
        case SubtractionConstant::b1:
            return helper(s, isospin, sub, evaluation, max_subs, 1., 1);
            
        default:
            throw std::domain_error{"Wrong input for subtraction constants."};
            
    }//switch
}


//-------------------------------------------------------------------------
//-------------------------------------------------------------------------
}
