#ifndef KERNEL_KHURI_HEADER
#define KERNEL_KHURI_HEADER

#include "cauchy.h"
#include "constants.h"
#include "curved_omnes.h"
#include "facilities.h"
#include "grid.h"
#include "helpers.h"
#include "mandelstam.h"
#include "omnes.h"
#include "phase_space.h"
#include "type_aliases.h"

#include "Eigen/Dense"

#include <functional>
#include <iterator>
#include <limits>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <type_traits>

/// @brief Solve KT equations via modified Gasser-Rusetsky method.

/// Provide means to solve KT equations for the scattering/decay involving
/// three pions with arbitrary mass and one particle with arbitrary mass and
/// quantum numbers
/// * I(JPC) = 0(1--);
/// * I(JPC) = 1(1-+).
/// The equations can be solved both iteratively and via direct matrix
/// inversion. The only facilities to be used directly are the class `Basis`
/// and the function `make_basis`, everything else may be consired as
/// implementation details.
namespace kernel {
using omnes::OmnesF;
using curved_omnes::CurvedOmnes;
using grid::Complex;
using grid::Grid;
using grid::Point;
using Matrix =
    Eigen::Matrix<Complex,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>;
using Vector = Eigen::VectorXcd;
using facilities::square;
using helpers::hits_threshold_m;
using type_aliases::CFunction;

/// The quantum numbers (JPC) of the decaying particle.
enum class Channel {
    one_minus_minus,
    one_minus_plus,
    two_plus_plus,
    eta
};

constexpr inline double angular_coefficient(Channel channel)
{
    switch (channel) {
        case Channel::one_minus_minus:
            return 1.5;
        case Channel::one_minus_plus:
            return -0.75;
        default:
            throw std::runtime_error{"unknown channel"};
    }
}

constexpr inline std::size_t index(std::size_t x_index, std::size_t z_index,
        std::size_t z_size)
    /// Convert from two-dimensional to one-dimensional indices.
{
    return x_index*z_size + z_index;
}

template<typename T>
inline double angular(const Grid<T>& g, std::size_t z_index)
    /// Compute angular contribution at given point of grid.
{
    return 1.0-square(g.z(z_index));
}

template<typename T>
Complex t_at(const Grid<T>& g, std::size_t x_index,
        std::size_t z_index, double pion_mass, double virtuality)
    /// Evaluate Mandelstam t at a given point of a grid.
{
    const auto p{g(x_index,z_index)};
    return mandelstam::t_photon_pion(p.x,p.z,pion_mass,virtuality);
}

template<typename F, typename T>
Vector sample_on_grid(const F& f, const Grid<T>& g, double pion_mass,
        double virtuality)
    /// Sample `f` at values of Mandelstam t on grid `g`.
{
    const std::size_t n_x{g.x_size()};
    const std::size_t n_z{g.z_size()};
    Vector result(n_x*n_z);
    for (std::size_t i{0}; i<n_x; ++i)
        for (std::size_t a{0}; a<n_z; ++a) {
            const Complex t{t_at(g,i,a,pion_mass,virtuality)};
            result(index(i,a,n_z)) = f(t);
        }
    return result;
}

inline double max_distance(const Vector& a, const Vector& b)
    /// Return the squared maximal entrywise difference of `a` and `b`.
{
    return (a-b).cwiseAbs2().maxCoeff();
}

template<typename T>
std::vector<Complex> generate_x_dependent(const OmnesF& o,
    const CFunction& pi_pi,
    const Grid<T>& g, double pion_mass, int subtractions)
    /// Generate the x_j dependent terms needed in the integration kernel.
{
    const std::size_t n_x{g.x_size()};
    std::vector<Complex> x_dependent(n_x);
    for (std::size_t j{0}; j<n_x; ++j) {
        const auto x{g.x(j)};
        x_dependent[j] = pi_pi(x)/o(x)*phase_space::sigma(pion_mass,x)
            /std::pow(x,subtractions);
    }
    return x_dependent;
}

template<typename T>
Matrix generate_kernel(const CurvedOmnes& o, const CFunction& pi_pi,
    const Grid<T>& g, double pion_mass, double virtuality, int subtractions,
    Channel channel)
    /// Compute the integration kernel.
{
    const std::size_t n_x{g.x_size()};
    const std::size_t n_z{g.z_size()};
    const std::size_t n{n_x*n_z};
    Matrix result(n,n);

    // x_j dependent terms
    const std::vector<Complex> x_dependent{
        generate_x_dependent(o.original(),pi_pi,g,pion_mass,subtractions)};

    // t(x_i,z_a) dependent terms
    std::vector<Complex> t(n);
    std::vector<Complex> t_dependent(n);
    for (std::size_t i{0}; i<n_x; ++i)
        for (std::size_t a{0}; a<n_z; ++a) {
            const std::size_t in{index(i,a,n_z)};
            t[in] = t_at(g,i,a,pion_mass,virtuality);
            t_dependent[in] = o(t[in])*std::pow(t[in],subtractions);
        }

    // create the matrix
    const double coeff{angular_coefficient(channel)/constants::pi()};
    for (std::size_t i{0}; i<n_x; ++i)
        for (std::size_t a{0}; a<n_z; ++a) {
            const std::size_t in{index(i,a,n_z)};
            const Complex t_term{t_dependent[in]};
            for (std::size_t j{0}; j<n_x; ++j) {
                const Complex x_term{x_dependent[j]};
                // `cauchy` is the only term that couples columns and rows.
                const Complex cauchy{g.x(j)-t[in]};
                for (std::size_t b{0}; b<n_z; ++b) {
                    const auto& point{g(j,b)};
                    const double weight{point.x_weight*point.z_weight};
                    result(in,index(j,b,n_z)) =
                        coeff*x_term*t_term*weight*angular(g,b)
                        *point.x_derivative/cauchy;
                }
            }
        }

    return result;
}

Vector iteration(const Matrix& kernel, const Vector& start, double accuracy,
        facilities::OnOffStream status=facilities::OnOffStream{});
    ///< @brief Solve KT equations iteratively.
    ///<
    ///< @param kernel the integration kernel
    ///< @param start the initial values for the basis function
    ///< @param accuracy the iteration terminates if this precision is reached
    ///< @param status in verbose mode, the number of the current iteration is
    ///< printed to the specified stream

Vector inverse(const Matrix& kernel, const Vector& start);
    ///< @brief Solve KT equations via matrix inversion.
    ///<
    ///< @param kernel the integration kernel
    ///< @param start the Omnes function times the subtraction polynomial

/// The different available solution methods.
enum class Method {
    iteration,
    inverse
};

class UnknownMethod : public std::exception {
public:
    const char* what() const noexcept override {return message.data();}
private:
    std::string message{"Unknown method."};
};

template<typename T>
std::vector<Vector> basis(const CurvedOmnes& o, const CFunction& pi_pi,
        int subtractions, const Grid<T>& g, double pion_mass,
        double virtuality, Method method=Method::inverse,
        std::optional<double> accuracy=std::nullopt,
        Channel channel=Channel::one_minus_minus)
    /// @brief Compute the set of basis vectors for a given KT problem.
    ///
    /// @param o the Omnes function
    /// @param pi_pi the pion pion scattering amplitude
    /// @param subtraction the number of subtractions
    /// @param g the grid on which the integrands of the KT equations
    /// are sampled
    /// @param pion_mass the pion mass
    /// @param virtuality the 'mass' squared of the I=0, J=1, P=C=-1
    /// particle. Might take on arbitrary values (i.e. 0 and negative
    /// values are alowed, too).
    /// @param method determine whether equations are solved iteratively
    /// or via direct matrix inversion
    /// @param accuracy allows to tune the accuracy of the solution if
    /// iteration is used.
{
    if (subtractions < 0) {
        throw std::runtime_error{"negative number of subtractions "};
    }
    Matrix kernel{
        generate_kernel(o,pi_pi,g,pion_mass,virtuality,subtractions, channel)};
    Vector omnes_start{sample_on_grid(o,g,pion_mass,virtuality)};

    std::vector<Vector> result;
    for (int i{0}; i<subtractions; ++i) {
        auto polynomial{[i](const Complex& s){return std::pow(s,i);}};
        Vector start{sample_on_grid(polynomial,g,pion_mass,virtuality)};
        start = start.cwiseProduct(omnes_start);
        switch (method) {
            case Method::iteration: {
                constexpr double default_value{1e-8};
                const double precision{accuracy ? *accuracy : default_value};
                result.push_back(iteration(kernel,start,precision));
                break; }
            case Method::inverse:
                result.push_back(inverse(kernel,start));
                break;
            default:
                throw UnknownMethod{};
        }
    }
    return result;
}

template<typename T>
/// The basis of the solution space to a KT equation.
class Basis {
public:
    Basis(const OmnesF& omn, const CFunction& pi_pi, int subtractions,
        const Grid<T>& g, double pion_mass, double virtuality,
        Method method=Method::inverse,
        std::optional<double> accuracy=std::nullopt,
        double minimal_distance=1e-4,
        Channel channel=Channel::one_minus_minus);
        ///< @param o the Omnes function
        ///< @param pi_pi the pion pion scattering amplitude
        ///< @param subtraction the number of subtractions
        ///< @param g the grid on which the integrands of the KT equations
        ///< are sampled
        ///< @param pion_mass the pion mass
        ///< @param virtuality the 'mass' squared of the I=0, J=1, P=C=-1
        ///< particle. Might take on arbitrary values (i.e. 0 and negative
        ///< values are allowed, too).
        ///< @param method determine whether equations are solved iteratively
        ///< or via direct matrix inversion
        ///< @param accuracy allows to tune the accuracy of the solution if
        ///< iteration is used.
    Complex operator()(std::size_t i, Complex s) const;
        ///< @brief Evaluate the basis function with subtraction polynomial
        ///< s^`i` at `s`.
    std::size_t size() const;
        ///< The number of basis functions (i.e., of subtractions).
    double virtuality() const noexcept { return virt; }
    double mass() const noexcept { return pion_mass; }

    const omnes::OmnesF& omnes() const
    {
        return curved_omn.original();
    }
private:
    gsl::Cquad integrate;

    CurvedOmnes curved_omn;
    std::vector<Vector> _basis;
    int subtractions;
    double pion_mass;
    double virt;
    double minimal_distance;
    Channel channel;

    Grid<T> grid;
    std::vector<cauchy::Interpolate> integrands;
};

template<typename T>
Basis<T> make_basis(const OmnesF& omn, const CFunction& pi_pi,
        int subtractions, const Grid<T>& g, double pion_mass,
        double virtuality, Method method=Method::iteration,
        std::optional<double> accuracy=std::nullopt,
        Channel channel=Channel::one_minus_minus)
    /// @brief Generate a basis of the solution space to a KT equation.
    ///
    /// The arguments match exactly those of the class `Basis`.
{
    return Basis<T>{omn,pi_pi,subtractions,g,pion_mass,virtuality,method,
        accuracy, channel};
}

template<typename T>
std::vector<Complex> discrete_basis_integrand(const OmnesF& o,
        const CFunction& pi_pi, const Vector& basis, const Grid<T>& g,
        double pion_mass)
    /// @brief Return the Mandelstam-s independent part of the integrand needed
    /// in the evaluation of a basis function.
{
    const std::size_t n_x{g.x_size()};
    const std::size_t n_z{g.z_size()};
    std::vector<Complex> result(n_x);

    for (std::size_t j{0}; j<n_x; ++j) {
        for (std::size_t b{0}; b<n_z; ++b)
            result[j] += angular(g,b) * basis(index(j,b,n_z)) * g(j,b).z_weight;
        const auto x{g.x(j)};
        result[j] *= pi_pi(x)*phase_space::sigma(pion_mass,x)/o(x);
    }
    return result;
}

template<typename T>
cauchy::Interpolate basis_integrand(const OmnesF& o,
        const CFunction& pi_pi, const Vector& basis, const Grid<T>& g,
        double pion_mass)
    /// @brief Return the interpolated Mandelstam-s independent part of the
    /// integrand needed in the evaluation of a basis function.
{
    const auto discrete_integrand{
        discrete_basis_integrand(o,pi_pi,basis,g,pion_mass)};
    return cauchy::Interpolate{g.x_parameter_values(),discrete_integrand,
        gsl::InterpolationMethod::linear};
}

template<typename T>
std::vector<cauchy::Interpolate> basis_integrands(const OmnesF& o,
        const CFunction& pi_pi, const std::vector<Vector>& basis,
        const Grid<T>& g, double pion_mass)
    /// @brief Return the interpolated Mandelstam-s independent parts of the
    /// integrands needed in the evaluation of an entire basis.
{
    const auto n{basis.size()};
    std::vector<cauchy::Interpolate> result;
    for (std::size_t i{0}; i<n; ++i)
        result.push_back(basis_integrand(o,pi_pi,basis[i],g,pion_mass));
    return result;
}

template<typename T>
Basis<T>::Basis(const OmnesF& omn, const CFunction& pi_pi,
        int subtractions, const Grid<T>& g, double pion_mass,
        double virtuality, Method method, std::optional<double> accuracy,
        double minimal_distance, Channel channel)
    :
    curved_omn{CurvedOmnes(omn, pi_pi, g)},
    _basis{basis(curved_omn,pi_pi,subtractions,g,pion_mass,virtuality,method,
                 accuracy, channel)},
    subtractions{subtractions},
    pion_mass{pion_mass},
    virt{virtuality},
    minimal_distance{minimal_distance},
    channel{channel},
    grid{g},
    integrands{basis_integrands(omn,pi_pi,_basis,grid,pion_mass)}
{
}

template<typename T, typename F>
Complex cut_prescription(Grid<T> grid, double lower, double upper, double s,
        F f, int subtractions, const gsl::Cquad& integrate)
    /// @brief Compute the dispersive integral with integrand `f` assuming that
    /// `s` hits the integration contour, i.e. via Cauchy principal value.
    ///
    /// Currently, this works only for linearly parametrised linear curves.
{
    const auto start{grid.curve_func(lower)};
    const auto end{grid.curve_func(upper)};
    const auto singularity{std::real((s-start) / (end-start))+lower};
    const auto fs{f(singularity)};
    const auto l{std::log((1.0-s/end) / (s/start - 1.0))};
    const auto sub{subtractions-1};
    auto h{
        [singularity
        ,fs
        ,sub
        ,s
        ,f = std::move(f)
        ,g = std::move(grid)
        ](double x)
        {
            const auto cx{g.curve_func(x)};
            return (f(x)/std::pow(cx,sub) - fs/std::pow(s,sub))
                /cx/(x-singularity);
        }};
    const auto result{
        std::get<0>(cauchy::c_integrate(h,lower,upper,integrate))};
    return std::pow(s,subtractions)*result
        + fs*(Complex{0.0,1.0}*constants::pi() + l);
}

template<typename T, typename F>
Complex ordinary_prescription(Grid<T> grid, double lower, double upper,
        const Complex& s, F f, int subtractions, const gsl::Cquad& integrate)
    /// @brief Compute the dispersive integral with integrand `f` assuming that
    /// `s` does not hit the integration contour.
{
    auto h{[subtractions,s,f = std::move(f),g = std::move(grid)](double x)
        {
            const auto cx{g.curve_func(x)};
            const auto dx{g.derivative_func(x)};
            return f(x)/std::pow(cx,subtractions)/(cx-s)*dx;
        }};
    const auto result{
        std::get<0>(cauchy::c_integrate(h,lower,upper,integrate))};
    return std::pow(s,subtractions)*result;
}

template<typename T>
constexpr bool tolerant_equal(T a, T b, T tolerance=1e-16)
    /// Check whether `a` and `b` are equal up to `tolerance`.
{
    return a-b<tolerance && b-a<tolerance;
}

template<typename Container>
void remove_element(Container& c, const typename Container::value_type& value)
    /// @brief Remove `value` from `c` if it is contained in it, otherwise
    /// leave `c` unchanged.
{
    const auto new_end{std::remove_if(c.begin(),c.end(),
            [value](const auto& v){return v==value;})};
    c.erase(new_end,c.cend());
}

template<typename T>
std::vector<std::pair<T,T>> segments(const std::vector<T>& points)
    /// Return all pairs of non-equal successive values.
{
    const auto s{points.size()};
    if (s<2)
        return {};

    std::vector<std::pair<T,T>> result;
    for (std::size_t i{0}; i<s-1; ++i) {
        const auto x{points[i]};
        const auto y{points[i+1]};
        if (!tolerant_equal(x,y))
            result.push_back(std::make_pair(x,y));
    }
    return result;
}

template<typename T>
std::vector<std::pair<T,T>> segments_without(const std::vector<T>& points,
        const std::pair<T,T>& value)
    /// Return all pairs of non-equal successive values excluding `value`.
{
    auto s{segments(points)};
    remove_element(s,value);
    return s;
}

template<typename T>
Complex Basis<T>::operator()(std::size_t i, Complex s) const
{
    if (hits_threshold_m(pion_mass,s,minimal_distance)) {
        const double shift{minimal_distance * 1.1};
        return ((*this)(i, s - shift) + (*this)(i, s + shift)) / 2.0;
    }
    const auto& integrand{integrands.at(i)};
    Complex dispersive_integral;
    if (const auto segment = grid.hits(s)) {
        const auto x0{grid.x_parameter_lower()};
        const auto x1{segment->first};
        const auto x2{segment->second};
        const auto x3{grid.x_parameter_upper()};
        const auto intervals{segments_without({x0,x1,x2,x3},*segment)};
        const auto sr{s.real()};
        dispersive_integral = cut_prescription(grid,x1,x2,sr,integrand,
                subtractions,integrate);
        for (const auto& i: intervals)
            dispersive_integral += ordinary_prescription(grid,i.first,i.second,
                    sr,integrand,subtractions,integrate);
    }
    else
        dispersive_integral = ordinary_prescription(
                grid,grid.x_parameter_lower(),grid.x_parameter_upper(),
                s,integrand,subtractions,integrate);

    return curved_omn(s)
        * (std::pow(s,i)
           + angular_coefficient(channel)/constants::pi()*dispersive_integral);
}

template<typename T>
std::size_t Basis<T>::size() const
{
    // The class invariant assures that `subtractions` is non-negative.
    using V = std::make_unsigned_t<decltype(subtractions)>;
    return static_cast<V>(subtractions);
}
} // kernel

#endif // KERNEL_KHURI_HEADER
