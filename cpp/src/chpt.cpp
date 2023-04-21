#include "chpt.h"

namespace chpt {
Complex sigma_fraction(double mass, const Complex& s)
{
    const Complex sig{sigma(mass,s)};
    return (1.0+sig) / (1.0-sig);
}

Complex log_sigma(double mass, const Complex& s)
{
    return log(sigma_fraction(mass,s));
}

Complex L_sigma(double mass, const Complex& s)
{
    const Complex sig{sigma(mass,s)};
    const Complex frac{1.0/sig};
    return square(frac) * (0.5*frac*log_sigma(mass,s) - 1.0);
}

Complex t4(double mass, Complex s, double pion_decay, double l_diff)
{
    const Complex ps = sigma(mass,s);
    const Complex pss = square(ps);
    const Complex lo = t2(mass,s,pion_decay);
    const Complex ls = L_sigma(mass,s);
    const Complex coeff = s*pss / (4608.0*pow(pi(),3)*pow(pion_decay,4));
    const Complex c_term = s * (l_diff + 1.0/3.0) - 7.5*square(mass);
    const Complex b_term = (pow(mass,4)*0.5/s * ((15.0 - 96.0*pss + 9.0*pow(ps,4))
                        * square(ls) - (146.0 - 50.0*pss) * ls + 41.0));
    const Complex imag = rho(mass,s)*square(lo);
    return coeff * (c_term - b_term) + Complex{0.0,1.0}*imag;
}
} // chpt
