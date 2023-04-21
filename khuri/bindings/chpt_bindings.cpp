#include "chpt.h"

#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

PYBIND11_MODULE(_khuri_chpt, m) {
    m.doc() = "The I=J=1 \\pi\\pi\\to\\pi\\pi ChPT-partial-wave amplitudes up"
              " to NLO in terms of the pion decay constant in the chiral limt.";

    m.def("lo", py::vectorize(chpt::t2),
        "The CHPT LO amplitude.",
        py::arg("mass"),
        py::arg("s"),
        py::arg("pion_decay"));

    m.def("nlo", py::vectorize(chpt::t4),
        "The CHPT NLO amplitude.",
        py::arg("mass"),
        py::arg("s"),
        py::arg("pion_decay"),
        py::arg("l_diff"));
}
