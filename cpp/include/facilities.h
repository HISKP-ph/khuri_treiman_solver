#ifndef FACILITIES_HEADER_GUARD
#define FACILITIES_HEADER_GUARD

#include <algorithm>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

/// Useful small facilities needed in many situations that are not specific
/// to Khuri Treiman equations/dispersion relations.
namespace facilities {

/// An interface for an output stream that can be switched on and off.

/// If 'on', `<<` has the expected result. If 'off', `<<` has no effect.
/// This is useful for instance if a function shall sometimes run in verbose
/// mode (writing status reports to the screen) and sometimes in quiet mode.
class OnOffStream {
public:
    OnOffStream() : status{true}, stream{std::cerr} {}
    OnOffStream(bool b, std::ostream& out) : status{b}, stream{out} {}

    void on() {status=true;}
    void off() {status=false;}

    template<class T>
    OnOffStream operator<<(const T& a)
    {
        if (status)
            stream<<a;
        return *this;
    }

private:
    bool status;
    std::ostream& stream;
};

template<class T>
constexpr T square(const T& x)
{
    return x*x;
}

template<class T>
constexpr T identity(const T& x)
{
    return x;
}

template<class Return, class Argument=Return>
class constant {
    const Return c;
public:
    constant(const Return& t) : c{t} {}
    const Return& operator()(const Argument& a) const {return c;}
};

template<class In, class Out>
constexpr Out identity(const In& x)
    /// @tparam In needs to be implicitly convertible into `Out`.
{
    return x;
}

template<class F, class G>
constexpr auto compose(F&& f, G&& g)
    /// Return the composition of `f` and `g`.

    /// Both `F` and `G` need to meet the requirements of Callable.
    /// Both need to take one argument only.
    /// The return type of `G` needs to be implicitly convertible to the
    /// argument type of `F`.
{
    return std::bind(f,std::bind(g,std::placeholders::_1));
}

template<class Number>
std::vector<Number> linspace(Number left, Number right, std::size_t size)
    /// @brief Return a `vector` containing `size` evenly spaced values in
    /// [`left`,`right`].

    /// If `left>right`, the values are in descending order.
    /// Since `left` and `right` are guaranteed to be contained in the `vector`,
    /// `size>=2` is required.
    ///
    /// @tparam Number should be a floating point type (e.g. double, float).
{
    if (size<2)
        throw std::invalid_argument{"linspace cannot create interval \
if size<2"};

    std::vector<Number> numbers(size);
    Number increment = (right-left) / (size-1);
    for (Number& x: numbers) {
        x = left;
        left += increment;
    }
    numbers.back() = right;
    return numbers;
}

inline std::ofstream open_write(const std::string& name, int precision=20)
    /// Open file named `name` for writing and set precision to `precision`.
{
    std::ofstream out{name};
    if (!out)
        throw std::runtime_error("could not open "+name);
    out<<std::scientific<<std::setprecision(precision);
    return out;
}

inline std::ifstream open_read(const std::string& name)
    /// Open file named `name` for reading.
{
    std::ifstream in{name};
    if (!in)
        throw std::runtime_error("could not open "+name);
    return in;
}
} // facilities


#endif // FACILITIES_HEADER_GUARD
