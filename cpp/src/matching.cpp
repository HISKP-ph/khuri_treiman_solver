#include "matching.h"


namespace match{

Matching::Matching(double mass_1, double mass_2, double mass_3, double mass_4, double epsilon, double validity)
:
mass_1{mass_1}, mass_2{mass_2}, mass_3{mass_3}, mass_4{mass_4},
epsilon{epsilon}, validity{validity},
threshold_12{pow(mass_1+mass_2,2)},  threshold_34{pow(mass_3+mass_4,2)},
pseudo_threshold_34{pow(mass_4-mass_3,2)}
{errors_matching();}


double Matching::parameter_simplification(const Interpolate& amplitude, double s, double i, double j, double k, double l, double m, double s1, double n, double p, double q){
    return (n*amplitude(s1)+i*amplitude(s) + j*epsilon*amplitude.derivative(s)
    + k*pow(epsilon,2.)*amplitude.derivative2(s) + p*amplitude(s+q*epsilon))/(l*pow(epsilon,m));
}

void Matching::errors_matching() {
    // Only allow for positive epsilon and validity to avoid unneccessary complications
    if (epsilon<0) {
        throw std::domain_error{"Negative 'epsilon' not allowed."};
    }
    if (validity<0.05) {
        throw std::domain_error{"Negative 'validity' not allowed. 'validity' needs to be greater than 0.05."};
    }
    if (mass_4<mass_1+mass_2+mass_3) {
        throw std::runtime_error{"Error: decay is kinematically not allowed."};
        }
    if (epsilon<validity) {
        throw std::runtime_error{"Error: matching point lies inside the range of validity."};
    }
    if (pseudo_threshold_34-validity<=threshold_12+validity){
        throw std::runtime_error{"Error: 'validity' is chosen to large. Overlapping range of validity."};
    }
    if (threshold_34-validity<=pseudo_threshold_34+validity){
        throw std::runtime_error{"Error: 'validity' is chosen to large. Overlapping range of validity."};
    }
}


//---------------------------------------------------------------------------------
// -- Matching the S-wave ---------------------------------------------------------
//---------------------------------------------------------------------------------
double Matching::s_match(double s, const Interpolate& amp){
    
        // Differentiate between every case for different s that can hit the poles within the dispersion integral:
        // s near scattering threshold threshold_12
    if (s<=threshold_12+validity && s>=threshold_12){
        double matching_point=threshold_12+epsilon; // match at a point in the vicinity of threshold_12
        // matching condition
        double a=parameter_simplification(amp, matching_point, 15,-12,4,8,1./2,0,0,0,0);
        double b=parameter_simplification(amp, matching_point,-5,8,-4,4,3./2,0,0,0,0);
        double c=parameter_simplification(amp, matching_point, 3,-4,4,8,5./2,0,0,0,0);

        // expansion valid in a small range around threshold_12
        return
        (a +b*(s-threshold_12) +c*pow(s-threshold_12,2.));
        }
    
        // s near threshold threshold_34
    else if (s<=threshold_34+validity && s>=threshold_34-validity){
            // due to square root: additionally distinguish for s below and above threshold_34,
            // one could also use s<threshold_34 and s>=threshold_34 without changing the result
        if (s<=threshold_34) {
            double matching_point=threshold_34-epsilon; // match at a point slightly below threshold_34
            // matching condition
            double a=parameter_simplification(amp, matching_point, 15,12,4,8,1./2,0,0,0,0);
            double b=parameter_simplification(amp, matching_point, -5,-8,-4,4,3./2,0,0,0,0);
            double c=parameter_simplification(amp, matching_point, 3,4,4,8,5./2,0,0,0,0);

            // expansion valid in a small range around threshold_34, for s<threshold_34
            return
            (a +b*(threshold_34-s) +c*pow(threshold_34-s,2.));

        }
        else if (s>threshold_34){
            double matching_point=threshold_34+epsilon; // match at a point slightly above threshold_34

            // matching condition
            double a=parameter_simplification(amp, matching_point, 15,-12,4,8,1./2,0,0,0,0);
            double b=parameter_simplification(amp, matching_point, 5,-8,4,4,3./2,0,0,0,0);
            double c=parameter_simplification(amp, matching_point, 3,-4,4,8,5./2,0,0,0,0);

            // expansion valid in a small range around threshold_34, for s>threshold_34
            return
            (a +b*(threshold_34-s) +c*pow(threshold_34-s,2.));
        }
    }
    return amp(s);
}

//---------------------------------------------------------------------------------
// -- Matching the P-wave ---------------------------------------------------------
//---------------------------------------------------------------------------------
double Matching::p_match(double s, const Interpolate& amp){
    
        // Differentiate between every case for different s that can hit the poles within the dispersion integral:
        // s near scattering threshold threshold_12
    if (s<=threshold_12+validity && s>=threshold_12){
        double matching_point=threshold_12+epsilon; // match at a point in the vicinity of threshold_12
        // matching condition
        double a=parameter_simplification(amp, matching_point, 35,-20,4,8,3./2,0,0,0,0);
        double b=parameter_simplification(amp, matching_point, -21,16,-4,4,5./2,0,0,0,0);
        double c=parameter_simplification(amp, matching_point, 15,-12,4,8,7./2,0,0,0,0);

        // expansion valid in a small range around threshold_12
        return (a +b*(s-threshold_12) +c*pow(s-threshold_12,2.));
    }

        // s near threshold threshold_34
    else if (s<=threshold_34+validity && s>=threshold_34-validity){
                // due to square root: additionally distinguish for s below and above threshold_34
                // one could also use s<threshold_34 and s>=threshold_34 without changing the result
        if (s<=threshold_34) {
            double matching_point=threshold_34-epsilon; // match at a point slightly below threshold_34
            // matching condition
            double a=parameter_simplification(amp, matching_point, 35,20,4,8,3./2,0,0,0,0);
            double b=parameter_simplification(amp, matching_point, -21,-16,-4,4,5./2,0,0,0,0);
            double c=parameter_simplification(amp, matching_point, 15,12,4,8,7./2,0,0,0,0);
            // expansion valid in a small range around threshold_34, for s<threshold_34
            return (a +b*(threshold_34-s) +c*pow(threshold_34-s,2.));
        }
        else if (s>threshold_34){
            double matching_point=threshold_34+epsilon; // match at a point slightly above threshold_34
            // matching condition
            double a=parameter_simplification(amp, matching_point, -35,20,-4,8,3./2,0,0,0,0);
            double b=parameter_simplification(amp, matching_point, -21,16,-4,4,5./2,0,0,0,0);
            double c=parameter_simplification(amp, matching_point, -15,12,-4,8,7./2,0,0,0,0);
            // expansion valid in a small range around threshold_34, for s>threshold_34
            // Note the minus sign!
            return -(a +b*(threshold_34-s) +c*pow(threshold_34-s,2.));
        }
    }
    return amp(s); // if s is not close to a pole reproduce the default amplitude
}

//---------------------------------------------------------------------------------
// -- Matching the D-wave ---------------------------------------------------------
//---------------------------------------------------------------------------------
double Matching::d_match(double s, const Interpolate& amp){
    
        // Differentiate between every case for different s that can hit the poles within the dispersion integral:
        // s near scattering threshold threshold_12
    if (s<=threshold_12+validity && s>=threshold_12){
        double matching_point=threshold_12+epsilon; // match at a point in the vicinity of threshold_12
        // matching condition
        double a=parameter_simplification(amp, matching_point, 63,-28,4,8,5./2,0,0,0,0);
        double b=parameter_simplification(amp, matching_point, -45,24,-4,4,7./2,0,0,0,0);
        double c=parameter_simplification(amp, matching_point, 35,-20,4,8,9./2,0,0,0,0);

        // expansion valid in a small range around threshold_12
        return (a +b*(s-threshold_12) +c*pow(s-threshold_12,2.));
    }

        // s near threshold threshold_34
    else if (s<=threshold_34+validity && s>=threshold_34-validity){
                // due to square root: additionally distinguish for s below and above threshold_34
                // one could also use s<threshold_34 and s>=threshold_34 without changing the result
        if (s<=threshold_34) {
            double matching_point=threshold_34-epsilon; // match at a point slightly below threshold_34
            // matching condition
            double a=parameter_simplification(amp, matching_point, 63,28,4,8,5./2,0,0,0,0);
            double b=parameter_simplification(amp, matching_point, -45,-24,-4,4,7./2,0,0,0,0);
            double c=parameter_simplification(amp, matching_point, 35,20,4,8,9./2,0,0,0,0);
            // expansion valid in a small range around threshold_34, for s<threshold_34
            return (a +b*(threshold_34-s) +c*pow(threshold_34-s,2.));
        }
        else if (s>threshold_34){
            double matching_point=threshold_34+epsilon; // match at a point slightly above threshold_34
            // matching condition
            double a=parameter_simplification(amp, matching_point, 63,-28,4,8,5./2,0,0,0,0);
            double b=parameter_simplification(amp, matching_point, 45,-24,4,4,7./2,0,0,0,0);
            double c=parameter_simplification(amp, matching_point, 35,-20,4,8,9./2,0,0,0,0);
            // expansion valid in a small range around threshold_34, for s>threshold_34
            // Note the minus sign!
            return -(a +b*(threshold_34-s) +c*pow(threshold_34-s,2.));
        }
    }
    return amp(s); // if s is not close to a pole reproduce the default amplitude
}
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------






//---------------------------------------------------------------------------------
// -- Functions needed for the dispersion integral --------------------------------
//---------------------------------------------------------------------------------

Complex Matching::s_wave_integrand_expansion(double s, const Interpolate& amp){
    // use this only in a small range around pseudo_threshold_34
    if (s<=pseudo_threshold_34) {
        double matching_point=pseudo_threshold_34-epsilon; // match at a point slightly below pseudo_threshold_34

        // matching conditions
        double b=parameter_simplification(amp, matching_point, 3,3,2,1,1./2,pseudo_threshold_34,-3,0,0);
        double c=parameter_simplification(amp, matching_point, -3,-4,-4,1,1.,pseudo_threshold_34,3,0,0);
        double d=parameter_simplification(amp, matching_point, 1,1,2,1,3./2,pseudo_threshold_34,-1,0,0);
        // expansion valid in a small range around pseudo_threshold_34, for s<pseudo_threshold_34
        return  b +c*sqrt(pseudo_threshold_34-s) + d*(pseudo_threshold_34-s);
    }
    else if (s>pseudo_threshold_34){
        double matching_point=pseudo_threshold_34+epsilon; // match at a point slightly above pseudo_threshold_34
        // matching condition
        double b=parameter_simplification(amp, matching_point, 3,-3,2,1,1./2,pseudo_threshold_34,-3,0,0);
        double c=parameter_simplification(amp, matching_point, 3,-4,4,1,1.,pseudo_threshold_34,-3,0,0);
        double d=parameter_simplification(amp, matching_point, -1,1,-2,1,3./2,pseudo_threshold_34,1,0,0);

        // expansion valid in a small range around pseudo_threshold_34, for s>pseudo_threshold_34
        return (-1.i*b + c *sqrt(std::complex<double>(pseudo_threshold_34-s))+ (-1.i)*d*(pseudo_threshold_34-s));
    }
//    throw std::domain_error{"'s_wave_integrand_expansion' can only be used within the range of matching."};
    return 0;
}



Complex Matching::p_wave_integrand_expansion(double s, const Interpolate& amp){
    // use this only in a small range around pseudo_threshold_34
    if (s<pseudo_threshold_34) {
            double matching_point=pseudo_threshold_34-epsilon; // match at a point slightly below pseudo_threshold_34
            // matching conditions
            double c=parameter_simplification(amp, matching_point, -8,-8,-4,1,3./2,pseudo_threshold_34,8,0,0);
            double d=parameter_simplification(amp, matching_point, 3,3,2,1,2,pseudo_threshold_34,-3,0,0);
            // expansion valid in a small range around pseudo_threshold_34, for s<pseudo_threshold_34
            return c + d*sqrt(pseudo_threshold_34-s);
        }
    else if (s>pseudo_threshold_34){
            double matching_point=pseudo_threshold_34+epsilon; // match at a point slightly above pseudo_threshold_34
            // matching condition
            double c=parameter_simplification(amp, matching_point, 8,-8,4,1,3./2,pseudo_threshold_34,-8,0,0);
            double d=parameter_simplification(amp, matching_point, 3,-3,2,1,2,pseudo_threshold_34,-3,0,0);
            // expansion valid in a small range around pseudo_threshold_34, for s>pseudo_threshold_34
            return -1.i* c + d*sqrt(std::complex<double>(pseudo_threshold_34-s));
        }
    else if (s==pseudo_threshold_34){
        double matching_point=pseudo_threshold_34;
            double c=parameter_simplification(amp, matching_point, -8,-8,-4,1,3./2,pseudo_threshold_34,8,0,0);
            return c;
    }
    return 0;
}

Complex Matching::d_wave_integrand_expansion(double s, const Interpolate& amp){
    // use this only in a small range around pseudo_threshold_34
    if (s<pseudo_threshold_34) {
            double matching_point=pseudo_threshold_34-epsilon; // match at a point slightly below pseudo_threshold_34
            // matching conditions
            double d=parameter_simplification(amp, matching_point, 28,24,18,7./2,5./2,pseudo_threshold_34,-27,-1,-3);
            double e=parameter_simplification(amp, matching_point, -56,-44,-40,28,3.,pseudo_threshold_34,53,3,-3);
            // expansion valid in a small range around pseudo_threshold_34, for s<pseudo_threshold_34
            return d + e*sqrt(pseudo_threshold_34-s);
        }
    else if (s>pseudo_threshold_34){
            double matching_point=pseudo_threshold_34+epsilon; // match at a point slightly above pseudo_threshold_34
            // matching condition
            double d=parameter_simplification(amp, matching_point, 28,-24,18,7./2,5./2,pseudo_threshold_34,-27,-1,3);
            double e=parameter_simplification(amp, matching_point, 56,-44,40,28,3.,pseudo_threshold_34,-53,-3,3);
            // expansion valid in a small range around pseudo_threshold_34, for s>pseudo_threshold_34
            return -1.i* d + e*sqrt(std::complex<double>(pseudo_threshold_34-s));
        }
    else if (s==pseudo_threshold_34){
        double matching_point=pseudo_threshold_34;
            double d=parameter_simplification(amp, matching_point, 28,24,18,7./2,5./2,pseudo_threshold_34,-27,-1,-3);
            return d;
    }
    return 0;
}

// this parameter is constant for s below and s above pseudo-threshold
double Matching::p_wave_parameter(double s, const Interpolate& amp){
    if (s<=pseudo_threshold_34) {
            // due to square root: additionally distinguish for s below and above pseudo_threshold_34
        double matching_point=pseudo_threshold_34-epsilon; // match at a point slightly below pseudo_threshold_34
            // matching conditions
        double b=parameter_simplification(amp, matching_point, 6,5,2,1,1.,pseudo_threshold_34,-6,0,0);
            // expansion valid in a small range around pseudo_threshold_34, for s<pseudo_threshold_34
        return b;
    }
    else{
        double matching_point=pseudo_threshold_34+epsilon; // match at a point slightly above pseudo_threshold_34
        // matching condition
        double b=parameter_simplification(amp, matching_point, -6,5,-2,1,1,pseudo_threshold_34,6,0,0);

        // expansion valid in a small range around pseudo_threshold_34, for s>pseudo_threshold_34
        return b;
    };
}

// this parameter is constant for s below and s above pseudo-threshold
std::tuple<double,double> Matching::d_wave_parameters(double s, const Interpolate& amp){
    if (s<=pseudo_threshold_34) {
            // due to square root: additionally distinguish for s below and above pseudo_threshold_34
        double matching_point=pseudo_threshold_34-epsilon; // match at a point slightly below pseudo_threshold_34
            // matching conditions
        double b=parameter_simplification(amp, matching_point, 112,80,32,28,1.,pseudo_threshold_34,-111,-1,-3);
        double c=parameter_simplification(amp, matching_point, -126,-114,-68,14,2.,pseudo_threshold_34,123,3,-3);
            // expansion valid in a small range around pseudo_threshold_34, for s<pseudo_threshold_34
        return {b,c};
    }
    else{
        double matching_point=pseudo_threshold_34+epsilon; // match at a point slightly above pseudo_threshold_34
        // matching condition
        double b=parameter_simplification(amp, matching_point, -112,80,-32,28,1.,pseudo_threshold_34,111,1,3);
        double c=parameter_simplification(amp, matching_point, -126,114,-68,14,2.,pseudo_threshold_34,123,3,3);

        // expansion valid in a small range around pseudo_threshold_34, for s>pseudo_threshold_34
        return {b,c};
    };
}
//---------------------------------------------------------------------------------
//---------------------------------------------------------------------------------
}
