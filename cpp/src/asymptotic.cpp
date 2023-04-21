#include "asymptotic.h"
using namespace asymptotic;

Asymptotic1::Asymptotic1(Function phasein, double matching_pointin, double limitin)
	: phase{phasein}, matching_point{matching_pointin}, limit{limitin}
{
	value = phase(matching_point);
	d_value = gsl::derivative(phase,matching_point,1e-4).first;
	param1 = pow(limit-value,2)/matching_point/d_value;
	param2 = (limit-value)/matching_point/d_value - 1.;
	if (param2 < -1.) {throw std::invalid_argument("The continuation of the phase"
		" to high energies exhibits a pole above the matching point.");}
}

Asymptotic1s::Asymptotic1s(Interpolate phasein, double matching_pointin, double limitin)
	: phase{phasein}, matching_point{matching_pointin}, limit{limitin}
{
	value = phase(matching_point);
	d_value = phase.derivative(matching_point);
	param1 = pow(limit-value,2)/matching_point/d_value;
	param2 = (limit-value)/matching_point/d_value - 1.;
	if (param2 < -1.) {throw std::invalid_argument("The continuation of the phase"
		" to high energies exhibits a pole above the matching point.");}
}

Asymptotic2::Asymptotic2(Function phasein, double lowerin, double upperin, double limitin)
	: phase{phasein}, lower{lowerin}, upper{upperin}, limit{limitin} {}

double Asymptotic2::connect(double s)
{
	double a = std::pow((s-lower),2);
	double b = 3.0*upper - 2.0*s - lower;
	double c = upper - lower;
	return a * b / std::pow(c,3);
}
