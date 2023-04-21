#include "cauchy.h"
#include "constants.h"
#include "gmock/gmock.h"
#include "gsl_interface.h"
#include "gtest/gtest.h"
#include "test_common.h"
#include <vector>

using cauchy::Complex;
using constants::pi;
using ::testing::ElementsAre;

Complex circle(double angle)
{
    const Complex i{0, 1.0};
    return std::exp(i * angle);
}

Complex circle_derivative(double angle)
{
    const Complex i{0, 1.0};
    return i * std::exp(i * angle);
}

Complex polynomial(const Complex& z)
{
    return std::pow(z, 2) + 3.0 * std::pow(z, 3);
}

TEST(Basic, Real)
{
    const std::vector<Complex> numbers{{1.0, 2.0}, {0.0, -3.1}, {-1.1, 20}};
    const auto real_parts{cauchy::real(numbers)};
    EXPECT_THAT(real_parts, ElementsAre(1.0, 0.0, -1.1));
}

TEST(Basic, Imag)
{
    const std::vector<Complex> numbers{{1.0, 2.0}, {0.0, -3.1}, {-1.1, 20}};
    const auto real_parts{cauchy::imag(numbers)};
    EXPECT_THAT(real_parts, ElementsAre(2.0, -3.1, 20));
}

TEST(Integrate, AlongRealAxis)
{
    const auto integrate{gsl::Cquad{}};
    const auto result{cauchy::c_integrate(circle, 0.0, 2.0 * pi(), integrate)};
    const auto value{std::get<0>(result)};
    constexpr double tolerance{1e-15};
    expect_near(value, {0.0, 0.0}, tolerance);
}

TEST(Integrate, AlongCircle)
{
    const auto integrate{gsl::Cquad{}};
    const auto result{cauchy::c_integrate(polynomial, circle, circle_derivative,
            0.0, pi(), integrate)};
    const auto value{std::get<0>(result)};
    constexpr double tolerance{1e-15};
    expect_near(value, {-2.0 / 3.0, 0.0}, tolerance);
}

TEST(Differentiation, Circle)
{
    constexpr double value{0.24};
    const auto computed{cauchy::derivative(circle, value, 1e-6)};
    expect_near(computed, circle_derivative(0.24), 1e-10);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
