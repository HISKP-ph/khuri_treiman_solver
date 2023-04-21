#ifndef PIECEWISE_H
#define PIECEWISE_H

#include "facilities.h"
#include "grid.h"
#include "mandelstam.h"

#include <algorithm>
#include <numeric>
#include <type_traits>
#include <vector>

/// @brief Facilities for the generation of a grid, whose associated curve in
/// the x-plane is piecewise linear.
namespace piecewise {
using facilities::square;
using grid::Complex;
using grid::Curve;

/// a piecewise linear path in the complex plane
class Piecewise : public Curve {
public:
    /// the different available parametrisations
    enum Para {
        linear,
        quadratic
    };

    class UnknownPara : public std::exception {
    public:
        const char* what() const noexcept {return message.data();}
    private:
        std::string message{"Invalid choice of parametrisation."};
    };

    template<typename T>
    static std::vector<Piecewise::Para> all_linear(T size)
    {
        return std::vector<Piecewise::Para>(size,linear);
    }

    Piecewise(const std::vector<Complex>& knots,
            const std::vector<Para>& parametrisations);
        ///< @param knots The points in the complex plane that are connected
        ///< with straigt lines. The curve starts at the first element of knots
        ///< end ends at the last.
        ///< @param parametrisations The parametrisations for the different
        ///< segments of the piecewise curve.
    Piecewise(const std::vector<Complex>& knots);
        ///< @param knots The points in the complex plane that are connected
        ///< with straigt lines. The curve starts at the first element of knots
        ///< end ends at the last.
    Complex curve_func(double x) const override;
        ///< @brief Evaluate the piecewise linear curve at x.
        ///<
        ///< x=k corresponds to the k-th point (k=0,1,...).
    Complex derivative_func(double x) const override;
    Segment hits(const Complex& s) const override;
    std::vector<double> boundaries() const override;
    double lower() const noexcept {return 0.0;}
        ///< Return the parameter value corresponding to the start of the curve.
    double upper() const noexcept {return pieces.size();}
        ///< Return the parameter value corresponding to the end of the curve.
    std::size_t piece_index(double x) const;
        ///< @brief Return the number of the segment corresponding to the
        ///< parameter value x.
private:
    using CC = std::pair<Complex,Complex>;
    std::vector<Para> parametrisations;
    std::vector<CC> pieces;
    std::vector<CC> adjacent;
};

/// linear curve along the real axis
class Real : public Piecewise {
public:
    Real(double threshold, double cut,
         const std::vector<double>& intermediate={});
        ///< @brief The curve starts at `threshold` end ends at `cut`.
        ///<
        ///< Optionally, `intermediate` might contain a list of points in
        ///< between.
};

/// curve for the decay as described in the paper by Gasser and Rusetsky
class VectorDecay : public Piecewise {
public:
    VectorDecay(double pion_mass, double virtuality, double cut);
};

/// @brief curve for arbitrary virtualities above the three-pion threshold and
/// arbitrary pion masses
class Adaptive : public Piecewise {
public:
    Adaptive(double pion_mass, double virtuality, double cut);
};
} // piecewise

#endif // PIECEWISE_H
