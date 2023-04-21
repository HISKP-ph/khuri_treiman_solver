#include "curved_omnes.h"
#include "piecewise.h"

#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/functional.h"
#include "pybind11/numpy.h"

#include <utility>

namespace py = pybind11;
using omnes::OmnesF;
using curved_omnes::CFunction;
using curved_omnes::CurvedOmnes;

template<typename C>
CurvedOmnes make_curved(omnes::OmnesF o, CFunction amplitude, C curve)
{
    return CurvedOmnes(std::move(o), std::move(amplitude), std::move(curve));
}

template<typename C>
void create_binding(py::module&m, const std::string& name)
{
    m.def(("make_curved_" + name).c_str(),
          make_curved<C>,
          "generate a curved Omnes function",
          py::arg("omnes"),
          py::arg("amplitude"),
          py::arg("curve"));
}

PYBIND11_MODULE(_khuri_curved_omnes, m) {
    m.doc() = "The Omnes function with cut along a somewhat general curve.";

    py::class_<CurvedOmnes>(m, "_CurvedOmnes")
        .def("__call__", py::vectorize(&CurvedOmnes::operator()),
             py::arg("mandelstam_s"))
        .def("original", &CurvedOmnes::original);

    create_binding<piecewise::Real>(m, "real");
    create_binding<piecewise::VectorDecay>(m, "vector_decay");
    create_binding<piecewise::Adaptive>(m, "adaptive");
}
