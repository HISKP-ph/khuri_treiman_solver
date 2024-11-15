#include "iterative_solution.h"

#include "pybind11/pybind11.h"
#include "pybind11/complex.h"
#include "pybind11/functional.h"
#include "pybind11/numpy.h"
#include "pybind11/stl.h"

#include <vector>

namespace py = pybind11;

class PyIteration : public iteration::Iteration{
    public:
        /* Inherit the constructors */
        using iteration::Iteration::Iteration;

        /* Trampoline (need one for each virtual function) */
        void solve_khuri_treiman() override {
            PYBIND11_OVERRIDE_PURE(
                void,                      /* Return type */
                iteration::Iteration,      /* Parent class */
                solve_khuri_treiman        /* Name of function in C++ (must match Python name) */
            );
        }

        void write_output(SubtractionConstant sub, int iteration, const std::string output_file) override {
            PYBIND11_OVERRIDE_PURE(
                void,                   /* Return type */
                iteration::Iteration,   /* Parent class */
                write_output,           /* Name of function in C++ (must match Python name) */
                sub,                    /* Argument(s) */
                iteration,
                output_file    
            );
        }

        Complex evaluate_omnes(double s, int isospin, Setting evaluation) override {
            PYBIND11_OVERRIDE_PURE(
                Complex,                   /* Return type */
                iteration::Iteration,   /* Parent class */
                evaluate_omnes,           /* Name of function in C++ (must match Python name) */
                s,                    /* Argument(s) */
                isospin,
                evaluation    
            );
        }
};

PYBIND11_MODULE(_khuri_pinocchio, m) {
    m.doc() = "KT solutions using the Pionocchio integration path. "
              "Currently implemented are eta(')-> 3Pi, eta'->eta 2pi, V -> 3Pi, X -> 3Pi, T -> 3Pi. ";

    const std::string init_docstring =
        "Parameters\n"
        "----------\n"
        "mass_X involved masses in the decay Process 4->1+2+3\n"
        "phase_I the complete path name of the phase shift with isospin I that shall be imported\n"
        "epsilon the value at which an expansion (around the pseudo-threshold) within the dispersion integral is matched to the exact value\n"
        "validity the range of validity within which an expansion (around the pseudo-threshold) within the dispersion integral is used\n"
        "epsilon_tilde the value at which an expansion (around the remaining singularities)  within the dispersion integral is matched to the exact value\n"
        "validity_tilde the range of validity within which an expansion (around the remaining singularities) around the pseudo-threshold within the dispersion integral is used\n"
        "cutoff numeric cutoff of the dispersion integral\n"
        "iterations number of iterations\n"
        "maximal number of subtractions\n"
        "output_XX path name to export the basis solutions for the subtraction constant XX";

    py::enum_<enums::SubtractionConstant>(m, "SubtractionConstant")
        .value("a0", enums::SubtractionConstant::a0)
        .value("b0", enums::SubtractionConstant::b0)
        .value("a1", enums::SubtractionConstant::a1)
        .value("g1", enums::SubtractionConstant::g1)
        .value("h1", enums::SubtractionConstant::h1)
        .value("b1", enums::SubtractionConstant::b1)
        .value("sub_a0", enums::SubtractionConstant::sub_a0)
        .value("sub_b0", enums::SubtractionConstant::sub_b0)
        .value("sub_c0", enums::SubtractionConstant::sub_c0)
        .value("sub_d0", enums::SubtractionConstant::sub_d0)
        .value("sub_a1", enums::SubtractionConstant::sub_a1)
        .value("sub_b1", enums::SubtractionConstant::sub_b1);
        

    py::enum_<enums::Setting>(m, "Setting")
        .value("above", enums::Setting::above)
        .value("below", enums::Setting::below)
        .value("egg", enums::Setting::egg);


    py::class_<iteration::Iteration, PyIteration>(m, "Iteration")
        .def(py::init<double,
                      double,
                      double,
                      double,
                      double,
                      double,
                      double,
                      double,
                      double,
                      int,
                      int>(),
             py::arg("mass_1"),
             py::arg("mass_2"),
             py::arg("mass_3"),
             py::arg("mass_4"),
             py::arg("epsilon"),
             py::arg("validity"),
             py::arg("epsilon_tilde"),
             py::arg("validity_tilde"),
             py::arg("cutoff")=1000,
             py::arg("iterations"),
             py::arg("max_subs"));


    py::class_<iteration::IterationV3Pi, iteration::Iteration>(m, "IterationV3Pi")
        .def(py::init<double,
                      double,
                      const std::string&,
                      double,
                      double,
                      double,
                      double,
                      double,
                      int,
                      int,
                      const std::string&,
                      const std::string&>(),
             init_docstring.c_str(),
             py::arg("mass_1"),
             py::arg("mass_4"),
             py::arg("phase"),
             py::arg("epsilon"),
             py::arg("validity"),
             py::arg("epsilon_tilde"),
             py::arg("validity_tilde"),
             py::arg("cutoff")=1000,
             py::arg("iterations"),
             py::arg("max_subs"),
             py::arg("output_a"),
             py::arg("output_b"))
        .def("solve_kt", &iteration::IterationV3Pi::solve_khuri_treiman)
        .def("__call__", py::vectorize(&iteration::IterationV3Pi::operator()),
             py::arg("s"),
             py::arg("isospin"),
             py::arg("sub"),
             py::arg("set"),
             py::arg("iteration"))
        .def("evaluate_omnes", py::vectorize(&iteration::IterationV3Pi::evaluate_omnes),
             py::arg("s"),
             py::arg("isospin"),
             py::arg("set"));


    py::class_<iteration::IterationX3Pi, iteration::Iteration>(m, "IterationX3Pi")
        .def(py::init<double,
                      double,
                      const std::string&,
                      double,
                      double,
                      double,
                      double,
                      double,
                      int,
                      int,
                      const std::string&,
                      const std::string&>(),
             init_docstring.c_str(),
             py::arg("mass_1"),
             py::arg("mass_4"),
             py::arg("phase"),
             py::arg("epsilon"),
             py::arg("validity"),
             py::arg("epsilon_tilde"),
             py::arg("validity_tilde"),
             py::arg("cutoff")=1000,
             py::arg("iterations"),
             py::arg("max_subs"),
             py::arg("output_a"),
             py::arg("output_b"))
        .def("solve_kt", &iteration::IterationX3Pi::solve_khuri_treiman)
        .def("__call__", py::vectorize(&iteration::IterationX3Pi::operator()),
             py::arg("s"),
             py::arg("isospin"),
             py::arg("sub"),
             py::arg("set"),
             py::arg("iteration"))
        .def("evaluate_omnes", py::vectorize(&iteration::IterationX3Pi::evaluate_omnes),
             py::arg("s"),
             py::arg("isospin"),
             py::arg("set"));


    py::class_<iteration::IterationT3Pi, iteration::Iteration>(m, "IterationT3Pi")
        .def(py::init<double,
                      double,
                      const std::string&,
                      double,
                      double,
                      double,
                      double,
                      double,
                      int,
                      int,
                      const std::string&,
                      const std::string&>(),
             init_docstring.c_str(),
             py::arg("mass_1"),
             py::arg("mass_4"),
             py::arg("phase"),
             py::arg("epsilon"),
             py::arg("validity"),
             py::arg("epsilon_tilde"),
             py::arg("validity_tilde"),
             py::arg("cutoff")=1000,
             py::arg("iterations"),
             py::arg("max_subs"),
             py::arg("output_a"),
             py::arg("output_b"))
        .def("solve_kt", &iteration::IterationT3Pi::solve_khuri_treiman)
        .def("__call__", py::vectorize(&iteration::IterationT3Pi::operator()),
             py::arg("s"),
             py::arg("isospin"),
             py::arg("sub"),
             py::arg("set"),
             py::arg("iteration"))
        .def("evaluate_omnes", py::vectorize(&iteration::IterationT3Pi::evaluate_omnes),
             py::arg("s"),
             py::arg("isospin"),
             py::arg("set"));


    py::class_<iteration::IterationEta3Pi, iteration::Iteration>(m, "IterationEta3Pi")
        .def(py::init<double,
                      double,
                      const std::string&,
                      const std::string&,
                      const std::string&,
                      double,
                      double,
                      double,
                      double,
                      double,
                      int,
                      int,
                      const std::string&,
                      const std::string&,
                      const std::string&,
                      const std::string&,
                      const std::string&>(),
             init_docstring.c_str(),
             py::arg("mass_1"),
             py::arg("mass_4"),
             py::arg("phase_0"),
             py::arg("phase_1"),
             py::arg("phase_2"),
             py::arg("epsilon"),
             py::arg("validity"),
             py::arg("epsilon_tilde"),
             py::arg("validity_tilde"),
             py::arg("cutoff")=1000,
             py::arg("iterations"),
             py::arg("max_subs"),
             py::arg("output_a0"),
             py::arg("output_b0"),
             py::arg("output_a1"),
             py::arg("output_g1"),
             py::arg("output_h1"))
        .def("solve_kt", &iteration::IterationEta3Pi::solve_khuri_treiman)
        .def("__call__", py::vectorize(&iteration::IterationEta3Pi::operator()),
             py::arg("s"),
             py::arg("isospin"),
             py::arg("sub"),
             py::arg("set"),
             py::arg("iteration"))
        .def("evaluate_omnes", py::vectorize(&iteration::IterationEta3Pi::evaluate_omnes),
             py::arg("s"),
             py::arg("isospin"),
             py::arg("set"));
}


