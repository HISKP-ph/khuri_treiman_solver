#include "gsl_interface.h"

#include "pybind11/pybind11.h"

namespace py = pybind11;

PYBIND11_MODULE(_khuri_gsl, m) {
    m.doc() = "Interface to the gsl library.";

    py::class_<gsl::Settings>(m, "Settings")
        .def(py::init<double, double, std::size_t>(),
             py::arg("absolute_precision") = 0.0,
             py::arg("relative_precision") = 1e-7,
             py::arg("space") = 1000);
}
