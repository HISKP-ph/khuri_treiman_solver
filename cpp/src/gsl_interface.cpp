#include "gsl_interface.h"
#include <stdexcept>
#include <tuple>

namespace gsl {
// -- Error handling ----------------------------------------------------------
const TurnOffGslErrorsHelper TurnOffGslErrors::t{};

void check(int status)
{
    std::string message{gsl_strerror(status)};
    switch (status) {
        case 0: // everything worked fine
            return;
        case GSL_ENOMEM:  // could not allocate enough space
            throw Allocation_error{message};
        case GSL_EDIVERGE: // integral is divergent or too slowly convergent
            throw Divergence_error{message};
        case GSL_EMAXITER: // maximum number of subdivisions in integration
                           // exceeded
            throw Subdivision_error{message};
        case GSL_EROUND: // roundoff error
            throw Roundoff_error{message};
        case GSL_ESING: // bad integrand behaviour (e.g. non-integrable
                        // singularity)
            throw Bad_integrand_error{message};
        case GSL_EDOM: // domain error
            throw Domain_error{message};
        default:
            throw Error{message};
    }
}


// -- general helpers ---------------------------------------------------------

double unwrap(double x, void* p)
    // Call *p with x, where *p should be a std::function<double(double)>.
    // `unwrap` provides the signature needed by the gsl integration routine.
{
    auto fp = static_cast<std::function<double(double)>*>(p);
    return (*fp)(x);
}

// Wrap a general function such that in can be used by the gsl routines.
// ATTENTION: instances of this class are passed via pointers to the gsl
// routines. In addition, they own a function.
// Hence to be safe one needs to ensure that the instance exists
// for the entire run time of the gsl routine, which is achieved best via
// first defining the instance and on a next line invoking the gsl-routine,
// instead of the dangerous `gsl_routine(&Wrapper{f}, ...)`.
class Wrapper {
public:
    explicit Wrapper(Function f) : function{std::move(f)}
    {
        // Use the `void*`, which points to parameters usually, to pass a
        // `Function` to the gsl routine. `unwrap` takes care of everything,
        // i.e., applies this function. This way, the gsl routine can be used
        // both with function objects and ordinary functions.
        gsl_func.function = unwrap;
        gsl_func.params = &function;
    }
    
    gsl_function* operator&()
    {
        return &gsl_func;
    }
private:
    Function function;
    gsl_function gsl_func;
};


// -- Integration: Gauss-Legendre  --------------------------------------------

GaussLegendre::GaussLegendre(std::size_t s)
    : _size{s}, table{gsl_integration_glfixed_table_alloc(s)}
{
}

GaussLegendre::GaussLegendre(const GaussLegendre& other)
    : _size{other.size()}, table{gsl_integration_glfixed_table_alloc(_size)}
{
}

GaussLegendre::GaussLegendre(GaussLegendre&& other) noexcept
    : _size{other.size()}, table{std::move(other.table)}
{
}

GaussLegendre& GaussLegendre::operator=(const GaussLegendre& other)
{
    resize(other.size());
    return *this;
}

GaussLegendre& GaussLegendre::operator=(GaussLegendre&& other) noexcept
{
    _size = other.size();
    table = std::move(other.table);
    return *this;
}

std::pair<double,double> GaussLegendre::point(double lower, double upper,
        std::size_t i) const
{
    if (i>=_size)
        throw std::out_of_range{"requested value exceeds number of knots"};
    double point{-1}, weight{-1}; // useless values to indicate failure
    gsl_integration_glfixed_point(lower,upper,i,&point,&weight,table.get());
    return std::make_pair(point,weight);
}

double GaussLegendre::operator()(const Function& f, double lower,
        double upper) const
{
    Wrapper wrapper{f};
    return gsl_integration_glfixed(&wrapper,lower,upper,table.get());
}

void GaussLegendre::resize(std::size_t s)
{
    _size = s;
    table = std::unique_ptr<gsl_integration_glfixed_table,GlfixedDeleter>{
            gsl_integration_glfixed_table_alloc(_size)};
}

std::size_t GaussLegendre::size() const noexcept
{
    return _size;
}

GaussLegendre::~GaussLegendre() noexcept
{
}


// -- Integration: adaptive routines ------------------------------------------

Qag::Qag(const Settings& set)
: absolute_precision{set.absolute_precision},
    relative_precision{set.relative_precision},
    limit{set.space},
    workspace{set.space}
{
}

Value Qag::operator()(const Function& f, double lower, double upper) const
{
    int sign{signed_interval(lower,upper) ? 1 : -1};

    Wrapper wrapper{f};

    double result{0.0};
    double error{0.0};

    bool lower_inf{std::isinf(lower)};
    bool upper_inf{std::isinf(upper)};

    if (lower_inf && upper_inf)
        call(gsl_integration_qagi,&wrapper,absolute_precision,
                relative_precision,limit,workspace.data(),&result,&error);
    else if (lower_inf)
        call(gsl_integration_qagil,&wrapper,upper,absolute_precision,
                relative_precision,limit,workspace.data(),&result,&error);
    else if (upper_inf)
        call(gsl_integration_qagiu,&wrapper,lower,absolute_precision,
                relative_precision,limit,workspace.data(),&result,&error);
    else
        call(gsl_integration_qags,&wrapper,lower,upper,absolute_precision,
                relative_precision,limit,workspace.data(),&result,&error);

    return Value{sign*result,error};
}

void Qag::reserve(std::size_t space)
{
    workspace.reserve(space);
    limit = space;
}

Cquad::Cquad(const Settings& set)
: absolute_precision{set.absolute_precision},
    relative_precision{set.relative_precision},
    workspace{set.space}
{
}

Value Cquad::operator()(const Function& f, double lower, double upper) const
{
    int sign{signed_interval(lower,upper) ? 1 : -1};

    bool lower_inf{std::isinf(lower)};
    bool upper_inf{std::isinf(upper)};

    // `gsl_integration_cquad` does not provide functions for the integration
    // of infinite intervals. Hence, the required change of variables is
    // performed explicitly.
    Function integrand{f};
    if (lower_inf && upper_inf) {
        integrand = [&f](double x){return (f((1-x)/x) + f((x-1)/x)) / (x*x);};
        lower = 0.0;
        upper = 1.0;
    }
    else if (lower_inf) {
        integrand = [&f,upper](double x){return f(upper+(x-1)/x) / (x*x);};
        lower = 0.0;
        upper = 1.0;
    }
    else if (upper_inf) {
        integrand = [&f,lower](double x){return f(lower+(1-x)/x) / (x*x);};
        lower = 0.0;
        upper = 1.0;
    }

    Wrapper wrapper{integrand};
    double result{0.0};
    double error{0.0};

    std::size_t evaluations{0}; // dummy variable for function call
    call(gsl_integration_cquad,&wrapper,lower,upper,absolute_precision,
            relative_precision,workspace.data(),&result,&error,&evaluations);
    return Value{sign*result,error};
}

void Cquad::reserve(std::size_t space)
{
    workspace.reserve(space);
}

// -- Interpolation -----------------------------------------------------------

Interpolate::Interpolate(const std::vector<double>& x,
        const std::vector<double>& y, InterpolationMethod m, bool tolerant)
    : x_data{x}, y_data{y}, method{m}, tolerant{tolerant}
{
    if (x_data.size()!=y_data.size())
        throw std::invalid_argument("x and y need to have the same size");
    if (x_data.size()<method.min_size())
        throw std::invalid_argument("not enough data points for the choosen \
interpolation method");
    spline = gsl_interp_alloc(method.get_method(),x_data.size());
    call(gsl_interp_init,spline,x_data.data(),y_data.data(),x.size());
}

Interpolate::Interpolate(const Interpolate& other)
    : x_data{other.x_data}, y_data{other.y_data}, method{other.method},
    tolerant{other.tolerant}
{
    spline = gsl_interp_alloc(method.get_method(),x_data.size());
    call(gsl_interp_init,spline,x_data.data(),y_data.data(),x_data.size());
}

Interpolate::Interpolate(Interpolate&& other)
    : x_data{std::move(other.x_data)}, y_data{std::move(other.y_data)},
    method{other.method}, tolerant{other.tolerant},
    acc{other.acc}, spline{other.spline}
{
    other.acc = nullptr;
    other.spline = nullptr;
}

Interpolate::Interpolate()
    : method{InterpolationMethod::linear}
{
    spline = gsl_interp_alloc(method.get_method(),0); 
}

Interpolate& Interpolate::operator=(const Interpolate& other)
{
    // free old data
    gsl_interp_free(spline);

    // copy
    x_data = other.x_data;
    y_data = other.y_data;
    method = other.method;
    tolerant = other.tolerant;
    spline = gsl_interp_alloc(method.get_method(),x_data.size());
    call(gsl_interp_init,spline,x_data.data(),y_data.data(),x_data.size());

    return *this;
}

Interpolate& Interpolate::operator=(Interpolate&& other)
{
    // free old data
    gsl_interp_free(spline);

    // copy
    x_data = std::move(other.x_data);
    y_data = std::move(other.y_data);
    method = other.method;
    spline = other.spline;
    tolerant = other.tolerant;

    // empty old object
    other.acc = nullptr;
    other.spline = nullptr;

    return *this;
}

double Interpolate::operator()(double x) const
{
    return evaluate(gsl_interp_eval_e, x);
}

double Interpolate::eval(double x) const
{
    return this->operator()(x);
}

double Interpolate::derivative(double x) const
{
    return evaluate(gsl_interp_eval_deriv_e, x);
}

double Interpolate::derivative2(double x) const
{
    return evaluate(gsl_interp_eval_deriv2_e, x);
}

Interpolate::~Interpolate() noexcept
{
    gsl_interp_free(spline);
    gsl_interp_accel_free(acc);
}

InterpolationMethod::InterpolationMethod(InterpolationMethod::Method m)
    : MethodTemplate(nullptr, gsl_interp_type_min_size)
{
    switch (m) {
        case InterpolationMethod::linear:
            method =  gsl_interp_linear;
            break;
        case InterpolationMethod::polynomial:
            method =  gsl_interp_polynomial;
            break;
        case InterpolationMethod::cubic:
            method =  gsl_interp_cspline;
            break;
        case InterpolationMethod::cubic_periodic:
            method =  gsl_interp_cspline_periodic;
            break;
        case InterpolationMethod::akima:
            method =  gsl_interp_akima;
            break;
        case InterpolationMethod::akima_periodic:
            method =  gsl_interp_akima_periodic;
            break;
        case InterpolationMethod::steffen:
            method =  gsl_interp_steffen;
            break;
    }
}

Interpolate sample(const Function& f, const Interval& i,
        InterpolationMethod m, bool tolerant)
{
    std::vector<double> y_values(i.size());
    std::transform(i.cbegin(),i.cend(),y_values.begin(),f);

    return Interpolate{i,y_values,m,tolerant};
}

// -- 2D Interpolation -----------------------------------------------------------

InterpolationMethod2D::InterpolationMethod2D(InterpolationMethod2D::Method2D m)
        : MethodTemplate(nullptr, gsl_interp2d_type_min_size)
{
    switch (m) {
        case InterpolationMethod2D::bilinear:
            method =  gsl_interp2d_bilinear;
            break;
        case InterpolationMethod2D::bicubic:
            method =  gsl_interp2d_bicubic;
            break;
    }
}

Interpolate2D::Interpolate2D(const std::vector<double>& x,
        const std::vector<double>& y, const std::vector<double>& z, InterpolationMethod2D m)
    : x_data{x}, y_data{y}, z_data{z}, method{m}
{
    if (x_data.size()*y_data.size()!=z_data.size())
        throw std::invalid_argument("x*y and z need to have the same size");
    if (x_data.size()*y_data.size()<method.min_size())
        throw std::invalid_argument("not enough data points for the choosen \
interpolation method");
    spline = gsl_spline2d_alloc(method.get_method(),x_data.size(),y_data.size());
    call(gsl_spline2d_init,spline,x_data.data(),y_data.data(),z_data.data(),x.size(),y.size());
}

Interpolate2D::Interpolate2D(const Interpolate2D& other)
    : x_data{other.x_data}, y_data{other.y_data}, z_data{other.z_data}, method{other.method}
{
    spline = gsl_spline2d_alloc(method.get_method(),x_data.size(),y_data.size());
    call(gsl_spline2d_init,spline,x_data.data(),y_data.data(),z_data.data(),x_data.size(),y_data.size());
}

Interpolate2D::Interpolate2D(Interpolate2D&& other)
    : x_data{std::move(other.x_data)}, y_data{std::move(other.y_data)}, z_data{std::move(other.z_data)},
    method{other.method}, accx{other.accx}, accy{other.accy}, spline{other.spline}
{
    other.accx = nullptr;
    other.accy = nullptr;
    other.spline = nullptr;
}

Interpolate2D& Interpolate2D::operator=(const Interpolate2D& other)
{
    // free old data
    gsl_spline2d_free(spline);

    // copy
    x_data = other.x_data;
    y_data = other.y_data;
    z_data = other.z_data;
    method = other.method;
    spline = gsl_spline2d_alloc(method.get_method(),x_data.size(),y_data.size());
    call(gsl_spline2d_init,spline,x_data.data(),y_data.data(),z_data.data(),x_data.size(),y_data.size());

    return *this;
}

Interpolate2D& Interpolate2D::operator=(Interpolate2D&& other)
{
    // free old data
    gsl_spline2d_free(spline);

    // copy
    x_data = std::move(other.x_data);
    y_data = std::move(other.y_data);
    z_data = std::move(other.z_data);
    method = other.method;
    spline = other.spline;

    // empty old object
    other.accx = nullptr;
    other.accy = nullptr;
    other.spline = nullptr;

    return *this;
}

double Interpolate2D::operator()(double x, double y) const
{
    return evaluate(gsl_spline2d_eval_e, x, y);
}

double Interpolate2D::eval(double x, double y) const
{
    return this->operator()(x,y);
}

double Interpolate2D::derivativex(double x, double y) const
{
    return evaluate(gsl_spline2d_eval_deriv_x_e, x, y);
}

double Interpolate2D::derivativey(double x, double y) const
{
    return evaluate(gsl_spline2d_eval_deriv_y_e, x, y);
}

double Interpolate2D::derivative2xx(double x, double y) const
{
    return evaluate(gsl_spline2d_eval_deriv_xx_e, x, y);
}

double Interpolate2D::derivative2yy(double x, double y) const
{
    return evaluate(gsl_spline2d_eval_deriv_yy_e, x, y);
}

double Interpolate2D::derivative2xy(double x, double y) const
{
    return evaluate(gsl_spline2d_eval_deriv_xy_e, x, y);
}

Interpolate2D::~Interpolate2D() noexcept
{
    gsl_spline2d_free(spline);
    gsl_interp_accel_free(accx);
    gsl_interp_accel_free(accy);
}

Interpolate2D sample_2d(const std::function<double(double, double)>& f, const Interval& x, const Interval& y,
        InterpolationMethod2D m)
{
    const std::size_t n_x = x.size();
    const std::size_t n_y = y.size();
    std::vector<double> z_values(n_x*n_y);
    for (std::size_t i=0; i<n_x; ++i)
    {
        for (std::size_t j=0; j<n_y; ++j)
        {
            z_values[j*n_x+i]=f(x[i],y[j]);
        }
    }

    return Interpolate2D{x,y,z_values,m};
}

// -- Differentiation  --------------------------------------------------------

Value derivative(const Function& f, double value, double step_size,
        DerivativeMethod method)
{
    Wrapper wrapper{f};
    double result{0.0};
    double abs_error{0.0};
    switch (method) {
        case DerivativeMethod::central:
            gsl_deriv_central(&wrapper, value, step_size, &result, &abs_error);
            break;
        case DerivativeMethod::forward:
            gsl_deriv_forward(&wrapper, value, step_size, &result, &abs_error);
            break;
        case DerivativeMethod::backward:
            gsl_deriv_backward(&wrapper, value, step_size, &result, &abs_error);
            break;
        default:
            throw std::runtime_error{"switch does not cover all cases"};
    }
    return Value{result, abs_error};
}
} // gsl
