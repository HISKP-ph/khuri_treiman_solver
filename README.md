
[![DOI](https://zenodo.org/badge/630925513.svg)](https://zenodo.org/badge/latestdoi/630925513)

# khuri

Khuri-Treiman equations for different decay processes. The python API uses the Gasser-Rusetsky method, while the c++ API uses the Pinocchio method.\
In the c++ API the following processes are currently implemented:\
eta(') -> 3 pi (also C-violating)\
eta' -> eta pi pi (also C-violating)\
V -> 3 pi (J^PC=1^--)\
X -> 3 pi (J^PC=1^-+)\
T -> 3 pi (J^PC=2^++)\
\
The code contains all tools for decays into 3 particles of 1 or 2 different masses. And numerical tools are implementet for processes up to D-waves. 

## Methods

[NOTICE](./NOTICE.md) contains the references that are used in this project.

## Citation

How to cite this project is noted in [CITATION](./CITATION.cff).

## License

A GNU General Public License v3.0 is provided in [LICENSE](./LICENSE).

## Documentation

To build the documentation of the python API, go to the directory [docs](./docs).
Subsequently enter:

    make html

Then open the file ``docs/_build/html/index.html`` in a web browser.

To build the documentation of the C++ API, got to the directory [cpp/docs](./cpp/docs).
Subsequently enter:

    doxygen Doxyfile

Then open the file ``cpp/docs/build/html/index.html`` in a web browser.
