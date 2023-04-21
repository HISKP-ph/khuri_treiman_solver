#include "iam.h"

#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

PYBIND11_MODULE(_khuri_iam, m) {
    m.doc() = "The I=J=1 \\pi\\pi->\\pi\\pi IAM-partial-wave amplitudes up to"
              " NNLO \n\nThese amplitudes are expressed in terms of the pion"
              " decay constant in the chiral limit F.";

    m.def("nlo", py::vectorize(iam::iam_nlo),
        "The IAM amplitude up to NLO on the first Riemann sheet.",
        py::arg("mass"),
        py::arg("s"),
        py::arg("pion_decay"),
        py::arg("l_diff"));
}
