#include "mandelstam.h"

#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/numpy.h"

namespace py = pybind11;

auto t_vector_decay(mandelstam::Complex s, double z, double mass,
        double virtuality)
{
    return mandelstam::t_photon_pion(s, z, mass, virtuality);
}

PYBIND11_MODULE(_khuri_mandelstam, m) {
    m.doc() = "Provide Mandelstam variables for a general four-particle"
              " process as well as simplified computation in several special"
              " cases.";

    m.def("t_vector_decay", py::vectorize(t_vector_decay),
        "Mandelstam t for scattering/decay of a vector and/into three scalars.",
        py::arg("mandelstam_s"),
        py::arg("cosine"),
        py::arg("mass"),
        py::arg("virtuality"));

    m.def("a_photon_pion", py::vectorize(mandelstam::a_photon_pion),
        "The a in (Mandelstam t) = a + z b for photon+pion->pion+pion.",
        py::arg("mandelstam_s"),
        py::arg("mass"),
        py::arg("virtuality"));

    m.def("b_photon_pion", py::vectorize(mandelstam::b_photon_pion),
        "The b in (Mandelstam t) = a + z b for photon+pion->pion+pion.",
        py::arg("mandelstam_s"),
        py::arg("mass"),
        py::arg("virtuality"),
        py::arg("schwarz")=false);
}
