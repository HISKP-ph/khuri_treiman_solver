#include "path.h"

#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/numpy.h"
#include "pybind11/functional.h"
#include "pybind11/stl.h"

namespace py = pybind11;

PYBIND11_MODULE(_khuri_path, m) {
    m.doc() = "The Pinocchio path.";

    py::class_<path::Interval>(m, "Interval")
        .def_readwrite("lower", &path::Interval::lower)
        .def_readwrite("upper", &path::Interval::upper);

    py::class_<path::Path>(m, "Path")
        .def("derivative", py::vectorize(&path::Path::derivative),
             py::arg("curve_parameter"),
             py::arg("step_size"));

    py::class_<path::PolarEgg>(m, "PolarEgg")
        .def(py::init<double, double>(),
             py::arg("decay_mass"),
             py::arg("pion_mass"))
        .def("__call__", py::vectorize(&path::PolarEgg::operator()),
             py::arg("phi"))
        .def("phi", &path::PolarEgg::phi,
             py::arg("mandelstam_variable"));

    py::class_<path::Pinocchio, path::Path>(m, "Pinocchio")
        .def(py::init<double, double, double, double>(),
             py::arg("decay_mass"),
             py::arg("pion_mass"),
             py::arg("cut")=std::numeric_limits<double>::infinity(),
             py::arg("epsilon")=0.0)
        .def("__call__", py::vectorize(&path::Pinocchio::operator()),
             py::arg("curve_parameter"))
        .def("curve_parameter", &path::Pinocchio::curve_parameter,
             py::arg("mandelstam_s"))
        .def("domain", &path::Pinocchio::domain)
        .def("parameter_boundaries", &path::Pinocchio::parameter_boundaries)
        .def("mandelstam_boundaries", &path::Pinocchio::mandelstam_boundaries);
}
