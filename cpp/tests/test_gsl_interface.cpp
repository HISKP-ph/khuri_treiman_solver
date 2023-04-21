#include "facilities.h"
#include "gsl_interface.h"
#include "gtest/gtest.h"
#include "asymptotic.h"
#include "iam.h"
#include "constants.h"
#include <cmath>
#include <vector>

using gsl::GaussLegendre;
using gsl::Interpolate;
using gsl::Interpolate2D;
using constants::pi;
using asymptotic::Asymptotic1;

void test_integration(const GaussLegendre& g)
{
    double calculated{g([](double x){return 1.0;},0.0,1.0)};
    constexpr double tolerance{1e-6};
    EXPECT_NEAR(calculated,1.0,tolerance);
}

TEST(GaussLegendre, CopyConstructor)
{
    constexpr std::size_t size{100};
    GaussLegendre g{size};
    GaussLegendre g2{g};
    EXPECT_EQ(g.size(),size);
    EXPECT_EQ(g2.size(),size);

    test_integration(g);
    test_integration(g2);
}

TEST(GaussLegendre, MoveConstructor)
{
    constexpr std::size_t size{100};
    GaussLegendre g{size};
    GaussLegendre g2{std::move(g)};
    EXPECT_EQ(g2.size(),size);
    test_integration(g2);
}

TEST(GaussLegendre, CopyAsignment)
{
    constexpr std::size_t size{100};
    GaussLegendre g{size};
    GaussLegendre g2{size*2};
    EXPECT_EQ(g.size(),size);
    EXPECT_EQ(g2.size(),size*2);
    g2 = g;
    EXPECT_EQ(g2.size(),g.size());

    test_integration(g);
    test_integration(g2);

}

TEST(GaussLegendre, MoveAsignment)
{
    constexpr std::size_t size{100};
    GaussLegendre g{size};
    GaussLegendre g2{size*2};
    EXPECT_EQ(g.size(),size);
    EXPECT_EQ(g2.size(),size*2);
    g2 = std::move(g);
    EXPECT_EQ(g2.size(),size);

    test_integration(g2);
}

TEST(GaussLegendre, Size)
{
    constexpr std::size_t size{100};
    GaussLegendre g{size};
    EXPECT_EQ(g.size(),size);
}

TEST(GaussLegendre, Integrate)
{
    constexpr std::size_t size{3};
    constexpr double tolerance{1e-2};
    GaussLegendre g{size};
    const auto f{[](double x){return 2.0*std::pow(x,5) - x*x + 3.5*x - 1.0;}};
    double calculated{g(f,-2.0,5.0)};
    double expected{5172.42};
    EXPECT_NEAR(calculated,expected,tolerance);

    calculated = g(f,5.0,-2.0);
    EXPECT_NEAR(calculated,-expected,tolerance);
}

TEST(GaussLegendre, Point)
{
    constexpr std::size_t size{2};
    GaussLegendre g{size};
    constexpr double lower{-1.0};
    constexpr double upper{1.0};
    std::vector<std::pair<double,double>> points_weights{
        {-1.0/std::sqrt(3.0), 1.0},
        {1.0/std::sqrt(3.0), 1.0},
    };
    for (std::size_t i{0}; i<size; ++i) {
        auto point{g.point(lower,upper,i)};
        auto expected{points_weights[i]};
        EXPECT_DOUBLE_EQ(point.first,expected.first);
        EXPECT_DOUBLE_EQ(point.second,expected.second);
    }
}

TEST(GaussLegendre, Resize)
{
    constexpr std::size_t size{3};
    constexpr std::size_t new_size{70};

    GaussLegendre g{size};
    g.resize(new_size);

    EXPECT_EQ(g.size(),new_size);

    test_integration(g);
}

TEST(Interpolate, Sample)
{
    const std::vector<double> knots{1,2,3,4,5};
    const auto f{[](double x){return 2.0*x;}};
    const Interpolate i{gsl::sample(f,knots,gsl::InterpolationMethod::linear)};
    constexpr double x{2.5};
    EXPECT_DOUBLE_EQ(f(x),i(x));
}

TEST(Interpolate, Derivatives)
{
    const std::vector<double> knots{facilities::linspace(4.0, 6.0, 50)};
    const auto f{[](double x){return 3.0*x*x;}};
    const Interpolate i{gsl::sample(f,knots,gsl::InterpolationMethod::cubic)};
    constexpr double x{5.0};
    constexpr double tolerance{1e-8};
    EXPECT_NEAR(30.0, i.derivative(x), tolerance);
    EXPECT_NEAR(6.0, i.derivative2(x), tolerance);
}

TEST(Interpolate, Throw)
{
    std::vector<double> knots{1};
    auto f{[](double x){return 2.0*x;}};
    ASSERT_THROW(gsl::sample(f,knots,gsl::InterpolationMethod::linear),
            std::invalid_argument);
}

TEST(Interpolate2D, Calc)
{
	const std::vector<double> x{0, 1};
	const std::vector<double> y{0, 1};
	const std::vector<double> z{0,1,1,0.5};
	Interpolate2D i{x,y,z,gsl::InterpolationMethod2D::bilinear};
	EXPECT_DOUBLE_EQ(0.1,i(0,0.1));
	EXPECT_DOUBLE_EQ(0.5,i(0.5,0));
}

TEST(Interpolate2D, Sample)
{
	const std::vector<double> x{facilities::linspace(0., 4., 50)};
	const std::vector<double> y{facilities::linspace(0., 4., 50)};
    const auto f{[](double x, double y){return 3.*x*x+5.*y*y;}};
	const Interpolate2D i{gsl::sample_2d(f,x,y,
                                         gsl::InterpolationMethod2D::bicubic)};
	constexpr double tolerance{1e-8};
    EXPECT_NEAR(8.,i(1.,1.),tolerance);
	EXPECT_NEAR(23.,i(1.,2.),tolerance);
}

TEST(Interpolate2D, Derivative)
{
	const std::vector<double> x{facilities::linspace(0., 4., 50)};
	const std::vector<double> y{facilities::linspace(0., 4., 50)};
    const auto f{[](double x, double y){return 3.*x*x+5.*y*y;}};
	const Interpolate2D i{gsl::sample_2d(f,x,y,
                                         gsl::InterpolationMethod2D::bicubic)};
    constexpr double tolerance{1e-8};
	EXPECT_NEAR(9.,i.derivativex(1.5,2.),tolerance);
	EXPECT_NEAR(20.,i.derivativey(1.5,2.),tolerance);
	EXPECT_NEAR(6.,i.derivative2xx(1.5,2.),tolerance);
	EXPECT_NEAR(10.,i.derivative2yy(1.5,2.),tolerance);
	EXPECT_NEAR(0.,i.derivative2xy(1.5,2.),tolerance);

}

TEST(Asymptotic1, smooth)
{
    gsl::Function f{[](double s){return std::arg(iam::iam_nlo(0.13957,s,0.09,6.0));}};
    constexpr double matching{1.};
    Asymptotic1 phase_cont(f,matching,pi());
    constexpr double tolerance{1e-4};
    EXPECT_DOUBLE_EQ(f(matching-0.5),phase_cont(matching-0.5));
    EXPECT_DOUBLE_EQ(f(matching),phase_cont(matching));
    EXPECT_NEAR(f(matching+0.1),phase_cont(matching+0.1),tolerance*100);
    EXPECT_NEAR(pi(),phase_cont(1e5),tolerance);
}

void test_differentiation(gsl::DerivativeMethod method)
{
    auto f{[](double x){return x*x;}};
    auto [value, error] = gsl::derivative(f, 3.0, 1e-8, method);
    constexpr double tolerance{1e-5};
    EXPECT_NEAR(value, 6.0, tolerance);
    constexpr double max_err_size{5e-5};
    EXPECT_NEAR(error, 0.0, max_err_size);
}

TEST(Differentiation, Central)
{
    test_differentiation(gsl::DerivativeMethod::central);
}

TEST(Differentiation, Forward)
{
    test_differentiation(gsl::DerivativeMethod::forward);
}

TEST(Differentiation, Backward)
{
    test_differentiation(gsl::DerivativeMethod::backward);
}

int main(int argc, char** argv)
{
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
