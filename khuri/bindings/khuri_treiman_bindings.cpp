#include "khuri_treiman.h"
#include "omnes.h"

#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/functional.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"

#include <vector>
#include <tuple>
#include <type_traits>

namespace py = pybind11;

using khuri_treiman::Basis;
using khuri_treiman::Channel;
using khuri_treiman::Complex;
using khuri_treiman::Curve;
using khuri_treiman::CFunction;
using khuri_treiman::Grid;
using khuri_treiman::Method;
using khuri_treiman::OmnesBasis;
using khuri_treiman::Piecewise;
using khuri_treiman::Point;
using khuri_treiman::PWave;
using khuri_treiman::PWaveKernel;

template<typename T>
void create_grid_binding(py::module& m, const std::string& type_name)
{
    using V = std::add_lvalue_reference_t<std::add_const_t<T>>;
    using G = Grid<T>;
    const std::string name = "Grid" + type_name;
    const std::string init_docstring =
        "Parameters\n"
        "----------\n"
        "t: The continuous curve in the x-plane.\n"
        "x_sizes: The number of knots along the (different segments of the)"
        " curve in the x-plane.\n"
        "z_size: The number of knots along the line in the z-plane.\n";
    py::class_<G, T>(m, name.c_str())
        .def(py::init<V, std::vector<std::size_t>, std::size_t>(),
             init_docstring.c_str())
        .def(py::init<V, double, std::size_t>())
        .def("__call__", py::vectorize(&G::operator()),
             py::arg("x_index"),
             py::arg("z_index"))
        .def("x_parameter_values", &G::x_parameter_values)
        .def("x", &G::x,
             py::arg("x_index"))
        .def("derivative", &G::derivative,
             py::arg("x_index"))
        .def("z", &G::z,
             py::arg("z_index"))
        .def("x_size", &G::x_size)
        .def("z_size", &G::z_size)
        .def("x_parameter_lower", &G::x_parameter_lower)
        .def("x_parameter_upper", &G::x_parameter_upper);
}

template<typename T>
void create_basis_binding(py::module& m, const std::string& type_name)
{
    using B = Basis<T>;
    using G = Grid<T>;
    const std::string name = "Basis" + type_name;
    const std::string init_docstring =
        "Parameters\n"
        "----------\n"
        "o: the Omnes function\n"
        "pi_pi: the pion pion scattering amplitude\n"
        "subtraction: the number of subtractions\n"
        "g: the grid on which the integrands of the KT equations are sampled\n"
        "pion_mass: the pion mass\n"
        "virtuality: the 'mass' squared of the I=0, J=1, P=C=-1"
        " particle. Might take on arbitrary values (i.e. 0 and negative"
        " values are alowed, too).\n"
        "method: determine whether equations are solved iteratively"
        " or via direct matrix inversion\n"
        "accuracy: allows to tune the accuracy of the solution if"
        " iteration is used.";
    const std::string call_docstring =
         "Evaluate the basis function with subtraction polynomial s^`i` at `s`";
    py::class_<B>(m, name.c_str())
        .def(py::init<const omnes::OmnesF&,
                      const CFunction&,
                      int,
                      const G&,
                      double,
                      double,
                      Method,
                      std::optional<double>,
                      double,
                      Channel>(),
             init_docstring.c_str(),
             py::arg("o"),
             py::arg("pi_pi"),
             py::arg("subtractions"),
             py::arg("g"),
             py::arg("pion_mass"),
             py::arg("virtuality"),
             py::arg("method")=Method::inverse,
             py::arg("accuracy")=std::nullopt,
             py::arg("minimal_distance")=1e-4,
             py::arg("channel")=Channel::one_minus_minus)
        .def("__call__", py::vectorize(&B::operator()),
             call_docstring.c_str(),
             py::arg("i"),
             py::arg("s"))
        .def("size", &B::size)
        .def("virtuality", &B::virtuality)
        .def("mass", &B::mass);
}

template<typename B>
void create_plain_partial_wave_binding(py::module& m,
                                       const std::string& type_name)
{
    {
        using P = PWave<B>;
        const std::string name = "PWave" + type_name;
        const std::string init_docstring =
            "Parameters\n"
            "----------\n"
            "basis: a set of basis functions";
        const std::string call_docstring =
             "Evaluate the partial wave with subtraction polynomial s^`i` at `s`";
        const std::string evaluate_docstring = "Evaluate the partial wave.";
        py::class_<P>(m, name.c_str())
            .def(py::init<const B&>(),
                 init_docstring.c_str(),
                 py::arg("basis"))
            .def("__call__", py::vectorize(&P::operator()),
                 call_docstring.c_str(),
                 py::arg("i"),
                 py::arg("s"))
            .def("evaluate", py::vectorize(&P::evaluate),
                 evaluate_docstring.c_str(),
                 py::arg("s"),
                 py::arg("subtraction_constants"));
    }
    {
        using P = PWaveKernel<B>;
        const std::string name = "PWaveKernel" + type_name;
        const std::string init_docstring =
            "Parameters\n"
            "----------\n"
            "basis: a set of basis functions\n"
            "sites: the sites between which the basis function is interpolated\n"
            "spectral_subtractions: the number of subtractions in the spectral representation\n";
        const std::string call_docstring =
             "Evaluate the partial wave with subtraction polynomial s^`i` at `s`";
        const std::string evaluate_docstring = "Evaluate the partial wave.";
        py::class_<P>(m, name.c_str())
            .def(py::init<const B&, const gsl::Interval&, int>(),
                 init_docstring.c_str(),
                 py::arg("basis"),
                 py::arg("sites"),
                 py::arg("spectral_subtractions"))
            .def("__call__", py::vectorize(&P::operator()),
                 call_docstring.c_str(),
                 py::arg("i"),
                 py::arg("s"))
            .def("evaluate", py::vectorize(&P::evaluate),
                 evaluate_docstring.c_str(),
                 py::arg("s"),
                 py::arg("subtraction_constants"))
            .def("hat", py::vectorize(&P::hat),
                 py::arg("i"),
                 py::arg("s")) ;
    }
}

template<typename T>
void create_partial_wave_binding(py::module& m, const std::string& type_name)
{
    using B = Basis<T>;
    create_plain_partial_wave_binding<B>(m, type_name);
}

template<typename T>
void create_bindings(py::module& m, const std::string& type_name,
        bool partial_wave=false)
{
    create_grid_binding<T>(m, type_name);
    create_basis_binding<T>(m, type_name);
    if (partial_wave) {
        create_partial_wave_binding<T>(m, type_name);
    }
}

void create_omnes_basis_bindings(py::module& m)
{
    py::class_<OmnesBasis>(m, "OmnesBasis")
        .def(py::init<omnes::OmnesF, int, double, double>())
        .def("__call__", py::vectorize(&OmnesBasis::operator()),
             "Evaluate the basis function with subtraction polynomial"
             " s^`i` at `s`",
             py::arg("i"),
             py::arg("s"))
        .def("size", &OmnesBasis::size)
        .def("virtuality", &OmnesBasis::virtuality)
        .def("mass", &OmnesBasis::mass);

    create_plain_partial_wave_binding<OmnesBasis>(m, "OmnesBasis");
}

PYBIND11_MODULE(_khuri_khuri_treiman, m) {
    m.doc() = "Khuri Treiman equations.";

    py::enum_<Channel>(m, "Channel", "The quantum numbers (JPC)"
                                      "of the decaying particle.")
        .value("one_minus_minus", Channel::one_minus_minus)
        .value("one_minus_plus", Channel::one_minus_plus)
        .value("two_plus_plus", Channel::two_plus_plus)
        .value("eta", Channel::eta);

    py::class_<Curve>(m, "Curve", "A curve in the complex plane.")
        .def("curve_func", py::vectorize(&Curve::curve_func),
             "Evaluate the curve at `x`.",
             py::arg("x"))
        .def("derivative_func", py::vectorize(&Curve::derivative_func),
             "Evaluate the derivative of the curve at `x`.",
             py::arg("x"))
        .def("hits", &Curve::hits,
             "Determine whether `s` hits the curve.\n\n"
             "If `s` lies on the curve, return the parameter values marking"
             " the beginning and the end of the segment that is hit.",
             py::arg("x"))
        .def("boundaries", &Curve::boundaries);

    py::class_<Piecewise, Curve> piecewise(m, "Piecewise",
                                           "a piecewise linear path in the"
                                           " complex plane");
    piecewise
        .def(py::init<const std::vector<Complex>&,
                      const std::vector<Piecewise::Para>&>())
        .def("lower", &Piecewise::lower,
             "Return the parameter value corresponding to the start of the"
             " curve.")
        .def("upper", &Piecewise::upper,
             "Return the parameter value corresponding to the end of the"
             " curve.")
        .def("piece_index", py::vectorize(&Piecewise::piece_index),
                "Return the number of the segment corresponding to the"
                " parameter value `x`.",
             py::arg("x"));

    py::enum_<Piecewise::Para>(piecewise, "Para")
        .value("linear", Piecewise::Para::linear)
        .value("quadratic", Piecewise::Para::quadratic)
        .export_values();

    py::enum_<Method>(m, "Method",
                      "The different available solution methods.")
        .value("iteration", Method::iteration)
        .value("inverse", Method::inverse);

    py::class_<khuri_treiman::Real, Piecewise>(m, "Real",
                                            "linear curve along the real axis")
        .def(py::init<double, double, const std::vector<double>&>(),
             py::arg("threshold"),
             py::arg("cut"),
             py::arg("intermediate")=std::vector<double>{});

    py::class_<khuri_treiman::VectorDecay, Piecewise>(m, "VectorDecay",
                                            "curve for the decay as described"
                                            " in the paper by Gasser and"
                                            " Rusetsky")
        .def(py::init<double, double, double>());

    py::class_<khuri_treiman::Adaptive, Piecewise>(m, "Adaptive",
                                            "curve for arbitrary virtualities"
                                            " above the three-pion threshold"
                                            " and arbitrary pion masses")
        .def(py::init<double, double, double>());

    py::class_<Point>(m, "Point",
                      "A point in the (x,z)-plane.")
        .def(py::init<Complex,double,Complex,double,double>())
        .def(py::init<std::tuple<Complex,double,Complex,double,double>>())
        .def_readwrite("x", &Point::x)
        .def_readwrite("x_weight", &Point::x_weight)
        .def_readwrite("x_derivative", &Point::x_derivative)
        .def_readwrite("z", &Point::z)
        .def_readwrite("z_weight", &Point::z_weight);

    m.def("partial_wave_kernel",
          py::vectorize(&wave::kernel),
          py::arg("mandelstam_s"),
          py::arg("x"),
          py::arg("pion_mass"),
          py::arg("virtuality"),
          py::arg("subtractions"),
          py::arg("asymptotic")=std::nullopt);

    m.def("legendre",
          py::vectorize(&wave::legendre));

    create_bindings<khuri_treiman::Real>(m, "Real", true);
    create_bindings<khuri_treiman::VectorDecay>(m, "VectorDecay");
    create_bindings<khuri_treiman::Adaptive>(m, "Adaptive");

    create_omnes_basis_bindings(m);
}
