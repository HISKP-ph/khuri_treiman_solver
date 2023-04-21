#ifndef ENUMS_H
#define ENUMS_H

/// Defines useful enumerations
namespace enums{

/// Change Processes in examples/KT.cpp
enum class Process{
    eta_3_pi,
    etap_3_pi,
    etap_eta_pi_pi,
    V_3_pi,
    X_3_pi,
    T_3_pi
};

/// @brief Introduce the subtraction constants.
/// Each subtraction constant will have its own set of basis solutions
enum class SubtractionConstant{
    // only for eta->3pi and eta'->3pi
    a0,
    b0,
    a1,
    g1,
    h1,
    
    //these ones for eta'-> eta pi pi
    sub_a0,
    sub_b0,
    sub_c0,
    sub_d0,
    sub_a1,
    sub_b1,

    // for V->3Pi
    b1
};

/// @brief Chose where to evaluate the
/// dispersion-integrals and Omnès functions
enum class Setting{
    // infinitesimal above or below the cut
    above,
    below,
    
    // only for eta->3pi: along polar egg
    egg,
    
    // these ones for eta'-> eta pi pi to avoid crossing the cut
    // and to avoid endpoint-singularities
    plus_upper_a,
    plus_upper_b,
    plus_lower_a,
    plus_lower_b,
    minus_upper_a,
    minus_upper_b,
    minus_lower_a,
    minus_lower_b,
    zero_upper_a,
    zero_upper_b,
    zero_lower_a,
    zero_lower_b
};

/// Endpoint of Pinocchio integration
enum class PinocchioEndpoints{
    lower,
    upper
};

/// Different angular averages for processes with two different masses in the final state
enum class TypeOfAverage{
    plus,
    minus,
    zero
};

/// Needed for angular_averages in processes with two different masses in the final state
enum class IntegrationRegion{
    a,
    b
};

}

#endif // ENUMS_H

