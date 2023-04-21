#ifndef TYPE_ALIASES_H
#define TYPE_ALIASES_H

#include <complex>
#include <functional>

///Define type aliases for easier access.
namespace type_aliases {
using Complex = std::complex<double>;
using Curve = std::function<Complex(double)>;
using CFunction = std::function<Complex(Complex)>;
} // type_aliases

#endif // TYPE_ALIASES_H
